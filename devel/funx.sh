function fullPath {
  echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
  return 0
}

function ask {
  read -s -n 1 -r -p "$1 " yn; echo
  if [[ $yn =~ ^[Yy]$ ]]
  then
    return 0
  else
    if [[ $2 -eq 1 ]]; then [[ $yn =~ ^[Nn]$ ]] && return 1 || return 0;
    else return 1; fi
  fi
}

function pause {
  prompt=$1
  if [[ "$prompt" == "" ]]; then prompt="Press any key to continue: "; fi
  read -s -n 1 -r -p "$prompt "; echo
}

function is_zipped {
  if [[ -r $1 ]]; then
    return $(file $1 | grep -q 'gzip')
  else return 1; fi
}