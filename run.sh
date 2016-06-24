#!/bin/bash
rm -f *.o;

make;

if [ "$#" -lt 3 ]; then
	./main input.png 32 16
elif [ "$#" -ge 3 ]; then
	./main "$1" "$2" "$3"
else
	./main
fi