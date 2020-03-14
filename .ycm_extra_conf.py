from os.path import join
#BOHMHOME = "/home/ahn/Dropbox/tdse/bohm/cpp/bohm"

from os import getenv
BOHMHOME = getenv("BOHMHOME")
BOHMINC = join(BOHMHOME, "include")

GSLHOME = getenv("GSLHOME")
GSLINC = join(GSLHOME, "include")

def Settings( **kwargs ):
  return {
    'flags': [ 
        '-x', 'c++', '-Wall', '-Wextra', '-Werror', 
        '-I', BOHMINC,
        '-I', GSLINC,
        ],
  }

