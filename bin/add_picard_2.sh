#!/bin/bash
sed -i '/^@/{s/ /\/2 /}' $@
