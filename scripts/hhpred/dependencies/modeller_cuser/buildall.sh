#!/bin/bash

declare -a files
files[1]=""
files[3]="52"
files[5]="54"
files[7]="56"

for i in "${files[@]}"
do
	gcc -shared -Wall -fPIC \
		`$1/bin/mod9.13 --cflags --libs` \
		`pkg-config --cflags glib-2.0` \
		-I/usr/include/python2.7 \
		cuser_form${i}.c cuser_form${i}_wrap.c \
		-o $1/modlib/modeller/_cuser_form${i}.so -lm -lmodeller
done;	
