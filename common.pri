# * Add path to sources folder to the include search paths.
#   Necessary for GCC to find the .h files of (in Designer) promoted widgets.
#   When you promote widgets you must specify the .h relative to the "Sources" folder.
# * Pull in additional include paths from the custom INCLUDEPATHS environment variable.
INCLUDEPATH += $${_PRO_FILE_PWD_}/../Sources $$(INCLUDEPATHS)

INCLUDEPATH += /usr/include/lensfun

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp=libiomp5

LIBS += -lraw
