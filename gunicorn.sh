#!/bin/sh
gunicorn app:app -w 2 --threads 2 -b 0.0.0.0:5005 --timeout 360000
