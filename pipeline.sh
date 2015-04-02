#!/bin/bash

function fastqgrep() {
  local f=$1
  shift
  sed -n "n;p;n;n" $f | grep $@
}

function fullPath {
  echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
  return 0
}


if [[ ! -d bin ]]; then
  mkdir -p bin
  for f in devel/*; do
    [[ -x $f ]] && ln -s $(fullPath $f) bin/$(basename $f)
  done
fi
export PATH=$(fullPath bin/):$PATH
unset fullPath
