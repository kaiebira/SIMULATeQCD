#!/bin/bash

scriptdir=$(cd $(dirname $0); pwd)/podman-build
topdir=$(cd $scriptdir/../; pwd)
cd $topdir

if [[ $(id -u) = 0 ]]; then
  echo "Please do not run $(basename $0) as root, it is designed to be run as a normal user account that has podman permissions."
  exit 1
fi

if [[ $# -eq 0 ]] ; then
    echo 'No arguments provided. Run with -h for help.'
    exit 0
fi

touch $topdir/podman-build/container-info-simulateqcd.txt 

# This function parses a YAML file and converts it to a Bash script.
# Parameters:
#  $1: The path to the YAML file to parse
#  $2: The prefix to use for Bash variable names
#
# The expected output of this function is a series of Bash variable
# declarations, based on the contents of the input YAML file. For example if
# your YAML file contains 'RHEL_VERSION: 9', this would produce 
# CONFIG_RHEL_VERSION="9"
# Original code taken from: https://stackoverflow.com/a/21189044/4427375
function parse_yaml {
    local prefix=$2
    # Define regular expressions to match various parts of YAML syntax
    local s='[[:space:]]*'
    local w='[a-zA-Z0-9_]*'
    local fs=$(echo @|tr @ '\034')
    # Use sed to replace various parts of the YAML with Bash syntax
    sed -ne "s|^\($s\):|\1|" \
         -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
         -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
    # Use awk to create Bash variables based on the YAML data
    awk -F$fs '{
        indent = length($1)/2;
        vname[indent] = $2;
        for (i in vname) {
            if (i > indent) {
                delete vname[i]
            }
        }
        if (length($3) > 0) {
            vn="";
            for (i=0; i<indent; i++) {
                vn=(vn)(vname[i])("_")
            }
            printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
        }
    }'
}

BUILD_ARG=
DETACH_ARG="-d"
DEBUG_ARG=""
DEBUG_ENV=""

opts=$(getopt \
  -n $(basename $0) \
  -o h \
  --longoptions "build" \
  --longoptions "detach" \
  --longoptions "nodetach" \
  --longoptions "verbose" \
  --longoptions "debug" \
  -- "$@")
if [[ $? -ne 0 ]]; then
  opts="-h"
  echo
fi

eval set -- "$opts"
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h)
      echo "Usage:  $(basename $0) (build|list)"
      echo
      echo "Args:"
      echo "  build       Build SIMULATeQCD from source"
      echo "  list        Lists all the possible build targets"
      exit 0
      ;;
    --detach)
      DETACH_ARG="-d"
      ;;
    --nodetach)
      echo "User requested we do not detach from containers."
      DETACH_ARG=""
      ;;
    --verbose)
      VERBOSE="-vvv"
      ;;
    --)
      shift
      break
      ;;
  esac
  shift
done

# Check if user inadvertently enabled podman as root
if systemctl is-enabled podman podman.socket; then
  echo "It looks like you are running podman as root. You will need to disable this before continuing with 'systemctl disable --now podman podman.socket'"
  exit
fi

# Check if podman is installed
if ! which podman; then
  echo "You need to install podman to continue."
  exit
fi

# Enable podman for user if it isn't already on with the remote socket present
# The test checks if the socket for podman exists. Use -e because it isn't a 
# regular file.
if ! test -e $(podman info --format '{{.Host.RemoteSocket.Path}}'); then
  echo "podman service for user not enabled. Running 'systemctl enable --now --user podman podman.socket'"
  systemctl enable --now --user podman podman.socket
fi

# Export docker host and set it to the filepath of podman's socket
export DOCKER_HOST=unix://$(podman info --format '{{.Host.RemoteSocket.Path}}')

# Make container user UID match calling user so that containers dont leave droppings we cant remove
USER_ID=$(id -u)
GROUP_ID=$(id -g)
USERNAME=simulateqcd
GROUPNAME=simulateqcd

echo "Group ID: ${GROUP_ID}"
echo "User ID: ${USER_ID}"

# Call parse_yaml to create Bash variables from the YAML file
eval "$(parse_yaml "$scriptdir/config.yml")"

# Check how many cores are available. Only run if CORES is not set
if [ -z "$CORES" ]; then
    CORES=$(grep -c ^processor /proc/cpuinfo)

    # Output the number of cores found
    echo "Found ${CORES} cores."
    
    # Prompt the user to press Enter to continue
    read -p "WARNING: You have not set the number of cores to use. WE WILL USE ALL OF THEM. Note that it is not unusual for your computer to run very slowly when doing this because 100% of its processing power will be dedicated to compilation. It is also not uncommon to run out of memory while doing this and for the build to crash when compiling everything if you do not have enough memory. Single compilation targets are usually ok. Press Enter to continue or Ctrl+C to abort." 
fi

# Read in the YAML file
echo "RHEL_VERSION=$RHEL_VERSION"
echo "CUDA_VERSION=$CUDA_VERSION"
RHEL_VERSION=$RHEL_VERSION
CUDA_VERSION=$CUDA_VERSION

# Create the build directory
mkdir -p ${OUTPUT_DIRECTORY}

case $1 in
  rm)
    podman rm simulateqcd
    ;;

  list)

      # This function captures all the possible build targets from the CMakeLists.txt file

      echo '####################################'
      echo '#          BUILD TARGETS           #'
      echo '####################################'

      # Define the regular expression to match lines with format
      # add_SIMULATeQCD_executable(memManTest src/testing/main_memManTest.cpp)
      # The capture group ([[:alnum:]_]+) matches one or more alphanumeric characters 
      # or underscores. The rest is a literal string.
      regex='add_SIMULATeQCD_executable\(([[:alnum:]_]+)'

      # Find all matching lines in the CMakeLists.txt file
      matches=$(grep -oE "$regex" $topdir/CMakeLists.txt)

      # Extract the capture group from each match and add it to the list
      list=()
      while read -r match; do
        if [[ $match =~ $regex ]]; then
          list+=("${BASH_REMATCH[1]}")
        fi
      done <<< "$matches"

      # Sort the resulting list alphabetically
      sorted_list=($(echo "${list[@]}" | tr ' ' '\n' | sort -u))

      # Print the resulting sorted list
      printf '%s\n' "${sorted_list[@]}"

    ;;

  build)

      # Run Arch-specific checks for all the special people running Arch out there
      if [[ "$(grep -oP '^ID_LIKE=\K\w+' /etc/*-release)" =~ (arch|manjaro) ]]; then

        echo "######################################################################"
        echo "#       WE NOTICED YOU ARE RUNNING ARCH. YOU ARE VERY COOL.          #"
        echo "######################################################################"                                    

        echo "Running all the Arch-specific tests for all the Arch-specific issues."

        # Check the value of kernel.unprivileged_userns_clone
        value=$(sysctl -n kernel.unprivileged_userns_clone)

        if [ $value -eq 1 ]; then
          echo "kernel.unprivileged_userns_clone is already set to 1"
        else
          echo "Error: kernel.unprivileged_userns_clone is set to $value instead of 1"
          echo "Please run 'sudo sysctl -w kernel.unprivileged_userns_clone=1' to update it."
          echo "Note: Using `sysctl -w` will only set the value for the current session. To make it permanent, add the following line to /etc/sysctl.conf: kernel.unprivileged_userns_clone=1"
          exit 1
        fi

        if ! command -v slirp4netns &> /dev/null
        then
            echo "Error: slirp4netns is not installed. For Arch-based systems this package is a requirement for rootless Podman. Please install the slirp4netns package with 'sudo pacman -Sy slirp4netns'."
            exit 1
        else
            echo "slirp4netns is installed."
        fi
      fi

      # Get the current user's username
      username=$(whoami)

      # Get the current user's UID and GID
      uid=$(id -u $username)
      gid=$(id -g $username)

      # Check if the current user already has subuid and subgid entries
      if grep -q "^$username:" /etc/subuid && grep -q "^$username:" /etc/subgid; then
        echo "subuid and subgid entries already exist for $username"
      else
        echo "Your user does not have subuid and subgid entries. This is required for rootless Podman to work properly."
        echo "You need to run 'sudo usermod --add-subuids 100000-165535 --add-subgids 100000-165535 $username && podman system migrate'"
        echo "If you are not an administrator you will need to contact your administrator to do this for you."
        echo "For details see see https://wiki.archlinux.org/title/Podman#Set_subuid_and_subgid"
        exit 1
      fi

      # Enable the podman service for the user if it isn't already on
      systemctl enable --now --user podman podman.socket

      # Delete any old containers still running on the server
      for container in simulateqcd
      do
        id=$(podman container ls -a --filter name=$container -q)
        if [[ -z $id ]]; then
          continue
        fi
        echo "Cleaning up old containers for $container: $id"
        echo -n "Stopping: "
        podman container stop $id
        echo -n "Removing: "
        podman container rm -v $id
      done

      # Make sure PROFILE is set and valid
      if [[ "$PROFILE" != "nvidia" && "$PROFILE" != "hip_nvidia" && "$PROFILE" != "hip_amd" ]]; then
          echo "Invalid value for PROFILE in configuration."
          echo "Allowed values are: nvidia, hip_nvidia, hip_amd."
          exit 1
      fi

      # Check if the RHEL_VERSION environment variable is set to "latest" and if it is get the latest version
      if [[ "${RHEL_VERSION}" == "latest" ]]; then
          # Set the URL
          url="https://developer.download.nvidia.com/compute/cuda/repos/"

          # Get the page content and filter for lines containing "rhel"
          content=$(curl -s "$url" | grep "rhel")

          # Parse the content for version numbers
          versions=$(echo "$content" | sed -nE 's/.*rhel([0-9]+).*/\1/p' | sort -unr)

          # Check which Rocky Linux version actually exists
          RHEL_VERSION=""
          available_versions=""
          for version in $versions; do
              # Check if Rocky Linux image exists for this version
              if podman manifest inspect docker.io/library/rockylinux:${version}-minimal >/dev/null 2>&1; then
                  if [[ -z "$RHEL_VERSION" ]]; then
                      RHEL_VERSION="$version"
                  fi
                  available_versions="${available_versions} ${version}"
              fi
          done

          # Exit with error if detection fails
          if [[ -z "$RHEL_VERSION" ]]; then
              echo "ERROR: Could not auto-detect latest RHEL version."
              echo "Available Rocky Linux versions: ${available_versions:-none found}"
              echo ""
              echo "Please edit podman-build/config.yml and set RHEL_VERSION"
              echo "Example:"
              echo "  RHEL_VERSION: 9"
              exit 1
          fi
          
          echo "Auto-detected latest RHEL version: $RHEL_VERSION"
      fi

      # Handle CUDA version for nvidia and hip_nvidia profiles
      if [[ "$PROFILE" == "nvidia" || "$PROFILE" == "hip_nvidia" ]]; then

          # Function to check if a Docker image exists on Docker Hub
          # Parameters:
          #   $1: The Docker image namespace/name (e.g., "nvidia/cuda")
          #   $2: The image tag to check
          # Returns:
          #   0 if the image exists, 1 otherwise
          function docker_image_exists() {
              local image="$1"
              local tag="$2"
              
              # Parse namespace and repository from image
              if [[ "$image" == *"/"* ]]; then
                  namespace="${image%/*}"
                  repo="${image#*/}"
              else
                  namespace="library"
                  repo="$image"
              fi
              
              # Query Docker Hub API v2 for the specific tag
              response=$(curl -s "https://hub.docker.com/v2/namespaces/${namespace}/repositories/${repo}/tags/${tag}")
              
              # Check if the response contains the tag (not an error message)
              if echo "$response" | grep -q '"name"'; then
                  return 0
              else
                  return 1
              fi
          }

          # Function to get all available CUDA tags for a specific Rocky Linux version
          # Parameters:
          #   $1: The RHEL/Rocky Linux version number (e.g., "8", "9")
          #   $2: Optional "silent" flag to suppress progress messages
          # Returns:
          #   A sorted list of CUDA versions available for the specified Rocky Linux version
          function get_available_cuda_rocky_tags() {
              local rhel_version="$1"
              local silent="$2"
              local page_size=100
              local next_url="https://hub.docker.com/v2/namespaces/nvidia/repositories/cuda/tags?page_size=${page_size}"
              local all_tags=""
              
              # Only show progress message if not in silent mode
              if [[ "$silent" != "silent" ]]; then
                  echo "Querying Docker Hub for available CUDA images..." >&2
              fi
              
              # Docker Hub API paginates results, so we need to iterate through all pages
              while [[ -n "$next_url" && "$next_url" != "null" ]]; do
                  response=$(curl -s "$next_url")
                  
                  # Extract tags matching pattern: *-devel-rockylinux{version}
                  tags=$(echo "$response" | grep -oP '"name":"[^"]*-devel-rockylinux'${rhel_version}'"' | sed 's/"name":"//;s/"//')
                  all_tags="${all_tags}${tags}"$'\n'
                  
                  # Get next page URL (Docker Hub uses \u0026 for & in JSON)
                  next_url=$(echo "$response" | grep -oP '"next":"[^"]*"' | sed 's/"next":"//;s/"//' | sed 's/\\u0026/\&/g')
              done
              
              # Extract CUDA version numbers, sort them, and remove duplicates
              echo "$all_tags" | grep -E '^[0-9]+\.[0-9]+\.[0-9]+-devel-rockylinux'${rhel_version}'$' | \
                  sed 's/-devel-rockylinux'${rhel_version}'$//' | \
                  sort -V | uniq
          }

          # Function to get the latest CUDA version for a specific RHEL version
          # Parameters:
          #   $1: The RHEL/Rocky Linux version number
          # Returns:
          #   The newest available CUDA version for the specified RHEL version
          function get_latest_cuda_for_rhel() {
              local rhel_version="$1"
              local cuda_versions=$(get_available_cuda_rocky_tags "$rhel_version" "silent")
              
              if [[ -z "$cuda_versions" ]]; then
                  return 1
              fi
              
              # Get the last (newest) version from the sorted list
              echo "$cuda_versions" | tail -n1
          }

          # Function to dynamically discover all Rocky Linux versions with CUDA support
          # This avoids hardcoding version numbers and automatically adapts to new releases
          # Returns:
          #   A reverse-sorted list of Rocky Linux version numbers (newest first)
          function get_all_rocky_versions_with_cuda() {
              local page_size=100
              local next_url="https://hub.docker.com/v2/namespaces/nvidia/repositories/cuda/tags?page_size=${page_size}"
              local all_versions=""
              
              # Iterate through all pages of Docker Hub results
              while [[ -n "$next_url" && "$next_url" != "null" ]]; do
                  response=$(curl -s "$next_url")
                  
                  # Extract Rocky Linux version numbers from all CUDA tags
                  # Pattern: extract the number after "rockylinux" in tags like "12.8.0-devel-rockylinux9"
                  versions=$(echo "$response" | grep -oP '"name":"[^"]*-devel-rockylinux[0-9]+' | \
                            grep -oP 'rockylinux[0-9]+' | \
                            sed 's/rockylinux//' | \
                            sort -n | uniq)
                  
                  for ver in $versions; do
                      all_versions="${all_versions} ${ver}"
                  done
                  
                  # Get next page URL
                  next_url=$(echo "$response" | grep -oP '"next":"[^"]*"' | sed 's/"next":"//;s/"//' | sed 's/\\u0026/\&/g')
              done
              
              # Return unique versions in reverse order (newest first)
              echo "$all_versions" | tr ' ' '\n' | sort -rn | uniq
          }

          # Function to validate and set RHEL version
          function validate_rhel_version() {
              local requested_version="$1"
              
              if [[ "$requested_version" == "latest" ]]; then
                  echo "Determining latest RHEL version with CUDA support..." >&2
                  # Get all available versions dynamically
                  local available_versions=$(get_all_rocky_versions_with_cuda)
                  
                  if [[ -z "$available_versions" ]]; then
                      echo "Error: Could not find any Rocky Linux CUDA images on Docker Hub" >&2
                      return 1
                  fi
                  
                  # Return the first (highest) version
                  echo "$available_versions" | head -n1
                  return 0
              else
                  # Validate specific version
                  echo "Validating RHEL version ${requested_version}..." >&2
                  if get_available_cuda_rocky_tags "$requested_version" "silent" | grep -q .; then
                      echo "$requested_version"
                      return 0
                  else
                      echo "Error: No CUDA images found for Rocky Linux ${requested_version}" >&2
                      echo "Checking available RHEL versions with CUDA support..." >&2
                      
                      # Get all available versions dynamically
                      local available_versions=$(get_all_rocky_versions_with_cuda)
                      
                      if [[ -n "$available_versions" ]]; then
                          echo "Available RHEL versions with CUDA support:" >&2
                          for version in $available_versions; do
                              echo "  - RHEL ${version}" >&2
                          done
                      else
                          echo "No Rocky Linux versions with CUDA support found on Docker Hub" >&2
                      fi
                      return 1
                  fi
              fi
          }

          # Validate and set RHEL version
          RHEL_VERSION=$(validate_rhel_version "$RHEL_VERSION")
          if [[ $? -ne 0 ]]; then
              exit 1
          fi
          echo "Using RHEL version: $RHEL_VERSION"

          # Handle CUDA version
          if [[ "$CUDA_VERSION" == "latest" ]]; then
              CUDA_VERSION=$(get_latest_cuda_for_rhel "$RHEL_VERSION")
              if [[ $? -ne 0 || -z "$CUDA_VERSION" ]]; then
                  echo "Error: Could not determine latest CUDA version for RHEL ${RHEL_VERSION}"
                  echo "Available CUDA versions for Rocky Linux ${RHEL_VERSION}:"
                  get_available_cuda_rocky_tags "$RHEL_VERSION" | sed 's/^/  - /'
                  exit 1
              fi
               echo "Auto-detected latest CUDA version: $CUDA_VERSION"
          else
              # Validate specific CUDA version
              if ! docker_image_exists "nvidia/cuda" "${CUDA_VERSION}-devel-rockylinux${RHEL_VERSION}"; then
                  echo "Error: CUDA image nvidia/cuda:${CUDA_VERSION}-devel-rockylinux${RHEL_VERSION} is not valid"
                  echo "Available CUDA versions for Rocky Linux ${RHEL_VERSION}:"
                  get_available_cuda_rocky_tags "$RHEL_VERSION" | sed 's/^/  - /'
                  exit 1
              fi
              echo "Validated CUDA version: $CUDA_VERSION"
          fi

          # Final validation and display
          cuda_image="nvidia/cuda:${CUDA_VERSION}-devel-rockylinux${RHEL_VERSION}"
          echo "Using Docker image: $cuda_image"
          
          # One final check to be absolutely sure
          if ! docker_image_exists "nvidia/cuda" "${CUDA_VERSION}-devel-rockylinux${RHEL_VERSION}"; then
              echo "Error: Final validation failed for ${cuda_image}"
              exit 1
          fi
          echo "✓ Image validated and available on Docker Hub"
      fi

      if [[ "$PROFILE" == "hip_nvidia" ]]; then
          BACKEND="hip_nvidia"
      fi

      if [[ "$PROFILE" == "hip_amd" ]]; then
          BACKEND="hip_amd"
      fi

      if [[ "$PROFILE" == "hip_nvidia" || "$PROFILE" == "hip_amd" ]]; then
          if [[ "$USE_GPU_P2P" != "OFF" ]]; then
              echo "Error: USE_GPU_P2P should be set to OFF when using the hip_nvidia or hip_amd profiles."
              exit 1
          fi
      fi

      if [[ "$PROFILE" == "nvidia" ]]; then

          export BUILDAH_ISOLATION=chroot

          echo "Running: podman build \
--tag simulateqcd/simulateqcd:latest \
--label name=simulateqcd \
--build-arg CORES=${CORES} \
--build-arg RHEL_VERSION=${RHEL_VERSION} \
--build-arg CUDA_VERSION=${CUDA_VERSION} \
--build-arg USERNAME=${USERNAME} \
--build-arg GROUPNAME=${GROUPNAME} \
--build-arg USER_ID=${USER_ID} \
--build-arg GROUP_ID=${GROUP_ID} \
--build-arg ARCHITECTURE=${ARCHITECTURE} \
--build-arg USE_GPU_AWARE_MPI=${USE_GPU_AWARE_MPI} \
--build-arg USE_GPU_P2P=${USE_GPU_P2P} \
--build-arg TARGET=${TARGET} \
--build-arg ADDITIONAL_CMAKE_OPTIONS=${ADDITIONAL_CMAKE_OPTIONS} \
--build-arg ADDITIONAL_MAKE_OPTIONS=${ADDITIONAL_MAKE_OPTIONS} \
-f $scriptdir/Dockerfile.cuda \
$topdir"

          podman build \
            --tag simulateqcd/simulateqcd:latest \
            --label name=simulateqcd \
            --build-arg CORES=${CORES} \
            --build-arg RHEL_VERSION=${RHEL_VERSION} \
            --build-arg CUDA_VERSION=${CUDA_VERSION} \
            --build-arg USERNAME=${USERNAME} \
            --build-arg GROUPNAME=${GROUPNAME} \
            --build-arg USER_ID=${USER_ID} \
            --build-arg GROUP_ID=${GROUP_ID} \
            --build-arg ARCHITECTURE=${ARCHITECTURE} \
            --build-arg USE_GPU_AWARE_MPI=${USE_GPU_AWARE_MPI} \
            --build-arg USE_GPU_P2P=${USE_GPU_P2P} \
            --build-arg TARGET=${TARGET} \
            --build-arg ADDITIONAL_CMAKE_OPTIONS=${ADDITIONAL_CMAKE_OPTIONS} \
            --build-arg ADDITIONAL_MAKE_OPTIONS=${ADDITIONAL_MAKE_OPTIONS} \
            -f $scriptdir/Dockerfile.cuda \
            $topdir

      elif [[ "$PROFILE" == "hip_nvidia" ]]; then

          echo "Running: podman build \
--tag simulateqcd/simulateqcd:latest \
--label name=simulateqcd \
--build-arg CORES=${CORES} \
--build-arg RHEL_VERSION=${RHEL_VERSION} \
--build-arg CUDA_VERSION=${CUDA_VERSION} \
--build-arg USERNAME=${USERNAME} \
--build-arg GROUPNAME=${GROUPNAME} \
--build-arg USER_ID=${USER_ID} \
--build-arg GROUP_ID=${GROUP_ID} \
--build-arg ARCHITECTURE=${ARCHITECTURE} \
--build-arg USE_GPU_AWARE_MPI=${USE_GPU_AWARE_MPI} \
--build-arg USE_GPU_P2P=${USE_GPU_P2P} \
--build-arg TARGET=${TARGET} \
--build-arg ADDITIONAL_CMAKE_OPTIONS=${ADDITIONAL_CMAKE_OPTIONS} \
--build-arg ADDITIONAL_MAKE_OPTIONS=${ADDITIONAL_MAKE_OPTIONS} \
--build-arg USE_HIP_AMD=${USE_HIP_AMD} \
--build-arg BACKEND=${BACKEND} \
-f $scriptdir/Dockerfile.nvidia_hip \
$topdir"

          podman build \
            --tag simulateqcd/simulateqcd:latest \
            --label name=simulateqcd \
            --build-arg CORES=${CORES} \
            --build-arg RHEL_VERSION=${RHEL_VERSION} \
            --build-arg CUDA_VERSION=${CUDA_VERSION} \
            --build-arg USERNAME=${USERNAME} \
            --build-arg GROUPNAME=${GROUPNAME} \
            --build-arg USER_ID=${USER_ID} \
            --build-arg GROUP_ID=${GROUP_ID} \
            --build-arg ARCHITECTURE=${ARCHITECTURE} \
            --build-arg USE_GPU_AWARE_MPI=${USE_GPU_AWARE_MPI} \
            --build-arg USE_GPU_P2P=${USE_GPU_P2P} \
            --build-arg TARGET=${TARGET} \
            --build-arg ADDITIONAL_CMAKE_OPTIONS=${ADDITIONAL_CMAKE_OPTIONS} \
            --build-arg ADDITIONAL_MAKE_OPTIONS=${ADDITIONAL_MAKE_OPTIONS} \
            --build-arg USE_HIP_AMD=${USE_HIP_AMD} \
            --build-arg BACKEND=${BACKEND} \
            -f $scriptdir/Dockerfile.nvidia_hip \
            $topdir

      else
        echo "Invalid PROFILE value: $PROFILE. PROFILE must be set to either 'nvidia' or 'hip_nvidia'"
        exit 1
      fi

      # Run the container with a no-op command to keep it running
      podman run --name simulateqcd -d simulateqcd/simulateqcd:latest tail -f /dev/null
      echo "Checking if the container is running..."

      if [[ -z $(podman inspect --format '{{.State.Running}}' simulateqcd | grep -i true) ]]; then
          echo  "It appears that the SIMULATeQCD (simulateqcd) container failed to deploy correctly. Look at the build output for details."
          exit 1
      fi

      podman cp simulateqcd:/build ${OUTPUT_DIRECTORY}

      # Remove dangling images (images that are not tagged)
      echo "Removing dangling images..."
      #podman rmi $(podman images -f "dangling=true" -q)

      echo "The build has finished and your binaries are located in ${OUTPUT_DIRECTORY}build."
      ;;  

  *)
    echo "Run $(basename $0) build to build SIMULATeQCD."
    exit 1
    ;;
esac
