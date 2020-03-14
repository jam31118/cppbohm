from os.path import join
from os import getenv

BOHMHOME = getenv("BOHMHOME")
BOHMINC = join(BOHMHOME, "include")

GSLHOME = getenv("GSLHOME")
GSLINC = join(GSLHOME, "include")

TDSEHOME = getenv("TDSEHOME")
TDSEINC = join(TDSEHOME, "include")

def Settings( **kwargs ):
  return {
    'flags': [ 
        '-x', 'c++', '-Wall', '-Wextra', '-Werror', 
        '-I', BOHMINC,
        '-I', GSLINC,
        '-I', TDSEINC
        ],
  }

