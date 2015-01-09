#!/bin/bash
sed -i '/^@/{s/ /\/1 /}' $@
