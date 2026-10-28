#!/usr/bin/env python3
#
# Script that uploads the docs to  https://sites.ecmwf.int/docs/ectrans/
#
# To install sites-toolkit:
#
#     pip install sites-toolkit --upgrade -index-url https://get.ecmwf.int/repository/pypi-all/simple

from sites.sdk.sites import Authenticator, Site

import argparse
from pathlib import Path

scripts = Path(__file__).parent.resolve()
ectrans_docs = scripts.parent.resolve()

parser = argparse.ArgumentParser(description='Publish documentation')
parser.add_argument('--token', type=str, default="")
parser.add_argument('--user', type=str, default="")
parser.add_argument('--password', type=str, default="")
parser.add_argument('--html', type=str, default=f"{ectrans_docs}/site")

args = parser.parse_args()

if args.user and args.password :
    print( "Authentication via user/password" )
    my_authenticator = Authenticator.from_credentials(username=args.user,password=args.password)
elif args.token :
    print( "Authentication via token" )
    my_authenticator = Authenticator.from_token(token=args.token)
else :
    print( "ERROR: no token, or user/password was provided" )
    from sys import exit
    exit(1)
    
# Create a Site instance
my_site = Site(space='docs', name='ectrans')

# Create a content_manager
content_manager = my_site.get_content_manager(authenticator=my_authenticator)

# Upload all the contents of a directory inside the content directory
content_manager.upload(local_path=args.html, recursive=True)
