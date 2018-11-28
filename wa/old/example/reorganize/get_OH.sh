#!/bin/bash

awk 'BEGIN {ORS=" "} NF==2 {print $2} /O/ {print $4, "\n"}' new.xyz > zO-t.txt
awk 'BEGIN {ORS=" "} NF==2 {print $2} /H/ {print $4, "\n"}' new.xyz > zH-t.txt
