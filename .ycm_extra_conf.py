from os.path import join
BOHMHOME = "/home/ahn/Dropbox/tdse/bohm/cpp/bohm"
BOHMINC = join(BOHMHOME, "include")

def Settings( **kwargs ):
  return {
    'flags': [ '-x', 'c++', '-Wall', '-Wextra', '-Werror', '-I', BOHMINC ],
  }

