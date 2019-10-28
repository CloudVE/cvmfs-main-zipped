from galaxy import eggs

import pkg_resources
pkg_resources.require( "bx-python" )

import logging, os, sys, time, sets, tempfile, shutil
import data
from galaxy import util
from galaxy.datatypes.sniff import *
from cgi import escape
import urllib
from bx.intervals.io import *
from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.binary import Binary


class Mrh( Binary ):
    """Class describing a MRH file"""
    file_ext = "mrh"

    def sniff( self, filename ):
        # MRH is a binary file type.
        # The first 8 bytes of any mrh file is '42054205'.
        try:
            header = open( filename, 'rb' ).read(8)
            if binascii.b2a_hex( header ) == binascii.hexlify( '42054205' ):
                return True
            return False
        except:
            return False

Binary.register_sniffable_binary_format("mrh", "mrh", Mrh)
