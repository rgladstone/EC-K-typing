#!/bin/bash
#use the first arguement to provide the gff file
 
#This is a fix for bakta input until panaroo is updated, first adding in the panaroo_id field
sed -i 's/;product=/;panaroo_id=;product=/g' "$1"
#and then the inference field
sed -i 's/;panaroo_id=;/;panaroo_id=;prepanaroo_inference=Unknown_inference;/g' "$1"
