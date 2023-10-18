#!/bin/bash

if [ $# -ne 1 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには1個の引数が必要です。" 1>&2
  exit 1
fi

docker run -it -d --name $1 -u `id -u`:`id -g` -v $(pwd):/usr/app/mpm_rust -w /usr/app/mpm_rust rust cargo run --release -p mlsmpm --example $1