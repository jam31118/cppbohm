from os.path import join, exists
from os import getenv

BOHMHOME = getenv("BOHMHOME")
BOHMINC = join(BOHMHOME, "include")
assert exists(BOHMINC)

GSLHOME = getenv("GSLHOME")
GSLINC = join(GSLHOME, "include")
assert exists(GSLINC)

TDSEHOME = getenv("TDSEHOME")
TDSEINC = join(TDSEHOME, "include")
assert exists(TDSEINC)

PARAMHOME = getenv("PARAMHOME")
PARAMINC = join(PARAMHOME, "include")
assert exists(PARAMINC)

def Settings( **kwargs ):
  return {
    'flags': [ 
        '-x', 'c++', '-Wall', '-Wextra', '-Werror', 
        '-I', BOHMINC,
        '-I', GSLINC,
        '-I', TDSEINC,
        '-I', PARAMINC,
        '-DFIELD',
        ],
  }

