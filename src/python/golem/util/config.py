# vim: ts=3:sw=3:expandtab
import sys
import datetime

import os
import os.path
import subprocess
import re
import tempfile
import shutil
import shlex
try:
   import ast
except ImportError:
   ast = None

import golem.util.path as gpath

class GolemConfigError(Exception):
   def __init__(self, value):
      Exception.__init__(self, value)
      self._value = value

   def __str__(self):
      return "GoSam Configuration Error: %s" % self._value

   def __repr__(self):
      return "%s(%r)" % ( self.__class__, self._value)

   def value(self):
      return self._value

class Property:
   """
   This class is a wrapper for properties including their
   type and a description.

   The use of a wrapper allows for future changes of the syntax
   without affecting the source code in many places.
   """
   def __init__(self, name, description, type=str, default=None,
         experimental=False, options=None, sep=None, hidden=False, seps=None):
      """
      Note:
      1) in the case of type=list sep encodes the
      delimiter character (';' or ',') if sep=None
      the comma is used.
      2) in the case of type="nested_list" seps is a list
      which encodes the delimiter characters if seps=None
      then the default seps=[';' ,','] is used.
      """
      self._name = name
      self._description = description
      self._type = type
      self._default = default
      self._experimental = experimental
      self._options = options
      self._sep = sep
      self._hidden = hidden
      self._seps = seps

   def _guess_correct(self, options, *given):
      result = []
      for word in given:
         min_lev = len(word)
         min_opt = None
         for option in options:
            dist = levenshtein(word, option)
            if dist < min_lev:
               min_opt = option
               min_lev = dist
         if min_opt is not None and min_lev <= 5:
            result.append(("If you meant '%s' instead of '%s', " +
               "please correct and rerun!") % (min_opt, word))
      return result

   def check(self, conf):
      result = []
   
      if self._type == int:
         if self._options is not None:
            value = conf.getProperty(str(self))

            if value is not None:
               try:
                  ivalue = int(value)
                  if value not in self._options:
                     result.append(
                        ("The value (%d) of option '%s'" +
                        " is not in the valid range.")
                        % (ivalue, self))
               except ValueError:
                  result.append(
                        ("The value (%r) of option '%s'" +
                        " is not an integer number.")
                        % (value, self))
   
      elif self._type == bool:
         bool_values = [
               "1", "true", ".true.", "t", ".t.", "yes", "y",
         		"0", "false", ".false.", "f", ".f.", "no", "n"]
         value = conf.getProperty(str(self))

         if value is not None:
            if value.lower() not in bool_values:
               result.append(
                     "The value (%r) of option '%s' is not a boolean literal."
                     % (value, self))
               result.extend(self._guess_correct(bool_values, value))
   
      elif self._type == list:
         if self._options is not None:
            odds = []
            default = self.getDefault()
            sep = self.getSep()
            sval = conf.getProperty(str(self))
            if sval is not None:
               if sep is None:
                  values = sval.split(',')
               else:
                  values = sval.split(sep)

               for value in values:
                  vls = value.lower().strip()
                  if vls == "":
                     continue
                  if vls not in self._options:
                     odds.append(value)
   
            if len(odds) > 0:
               result.append("The option '%s' contains unexpected values: %s"
                     % (self, ", ".join(map(repr, odds))))
               result.extend(self._guess_correct(self._options, *odds))
      else:
         if self._options is not None:
            value = conf.getProperty(str(self))
            if value is not None:
               if value.lower() not in self._options:
                  result.append("Unexpected value (%r) for option '%s'."
                        % (value, self))
                  result.extend(self._guess_correct(self._options, value))

      return result

   def isExperimental(self):
      return self._experimental

   def isHidden(self):
      return self._hidden


   def getName(self):
      return self._name

   def getDescription(self):
      return self._description

   def getType(self):
      return self._type

   def getDefault(self):
      return self._default

   def getSep(self):
      return self._sep

   def getSeps(self):
      return self._seps

   def __str__(self):
      return self.getName()

   def __repr__(self):
      return "Property(%r, %r, %r, %r)" % (
            self._name, self._description,
            self._type, self._default)

class Properties:
   """
   This class provides a simplistic replacement for
   the Java class java.util.Properties.

   It is simplistic with respect to the limited number
   of escape sequences which are implemented:
     '\\', '\n', '\r', '\f', '\t'.
   As in the Java implementation, a single backslash
   in front of any other character is removed.
   """

   def __init__(self, defaults=None, **values):
      self._defaults = defaults
      self._map = {}
      for key, value in values.items():
         self.setProperty(key, str(value))

      self._decode = False

      self.cache = {}

   def decode(self):
      self._decode = True

   def nodecode(self):
      self._decode = False

   def getProperty(self, key, default=None):
      result = self._getProperty(key)
      if result is None:
         return default
      else:
         return result

   def _getProperty(self, key):
      if isinstance(key, Property):
         type = key.getType()
         name = str(key)
         default = key.getDefault()
         sep = key.getSep()
         seps = key.getSeps()
         if(type == int):
            return self.getIntegerProperty(name, default)
         elif(type == bool):
            return self.getBooleanProperty(name, default)
         elif(type == list):
            if sep is None:
               return self.getListProperty(name, default, ',')
            else:
               return self.getListProperty(name, default, sep)
         elif(type == "nested_list"):
            if seps is None:
               return self.getNestedListProperty(name, default, [';',','])
            else:
               return self.getNestedListProperty(name, default, seps)
         else:
            return self.getProperty(name, default)
      else:
         if key in self._map:
            result = self._map[key]
         elif self._defaults is not None:
            if key in self._defaults:
               result = self._defaults[key]
            else:
               return None
         else:
            return None

         if self._decode and isinstance(result, str):
            # won't work in python3:
            # return result.decode("string_escape")
            # This is not 100% correct but should work reasonably well:
            if ast is not None:
               return ast.literal_eval("'"+result+"'")
            else:
               return result.decode("string_escape")
         else:
            return result


   def getListProperty(self, key, default=None, delimiter=','):
      name = str(key)
      if name in self:
         value = self[name].split(delimiter)
         return list(map(lambda x: x.strip(), value))
      else:
         if default:
            return default.split(delimiter)
         else:
            return []

   def getNestedListProperty(self, key, default=None, delimiters=[';',',']):

      def split_recursive(string, delimiters, index=0):
         result = []
         try:
            current_delimiter = delimiters[index]
         except IndexError:
            return string.strip()
         for item in string.split(current_delimiter):
            result.append(split_recursive(item, delimiters, index+1))
         return result

      name = str(key)
      if name in self:
         return split_recursive(self[name], delimiters)
      else:
         if default:
            return split_recursive(default, delimiters)
         else:
            return []

   def getBooleanProperty(self, key, default=False):
      true_values = ["1", "true", ".true.", "t", ".t.", "yes", "y"]
      name = str(key)
      if name in self:
         value = self.getProperty(name, default=default).strip().lower()
         return value in true_values
      else:
         return default

   def getIntegerProperty(self, key, default=None):
      name = str(key)
      if name in self:
         value = self.getProperty(name, default=default).strip()
         try:
            return int(value)
         except ValueError as ex:
            raise GolemConfigError(
               "Property '%s' does not contain an integer value ('%s')." %
               (key, value))
      else:
         return default

   def copyProperties(self, other, *keys):
      for key in keys:
         self.setProperty(key, other.getProperty(key))

   def copy(self):
      result = Properties()
      for key in self.propertyNames():
         result.setProperty(key, self.getProperty(key))
      return result

   def setProperty(self, key, value):
      name = str(key)
      if name.startswith("+"):
         self.setProperty(name[1:], value)
      if value.__class__ == list:
         self._map[name] = ",".join(map(str, value))
      else:
         self._map[name] = str(value)

   def __contains__(self, key):
      name = str(key)
      if name in self._map:
         return True
      elif self._defaults is not None:
         return name in self._defaults
      else:
         return False
         
   def propertyNames(self):
      for key in self._map.keys():
         yield key
      if self._defaults is not None:
         for key in self._defaults:
            if key not in self._map:
               yield key

   def list(self, stream=sys.stdout):
      for key in self:
         stream.write("%s=%s\n" % (escape(key, True), escape(self[key])))

   def store(self, stream, comments=None, properties=None,
         info=[]):

      def format_comment(prop):
         if prop.getType() == str:
            stype = "text"
         elif prop.getType() == int:
            stype = "integer number"
         elif prop.getType() == bool:
            stype = "true/false"
         elif prop.getType() == list:
            stype = "comma separated list"
         elif prop.getType() == "nested_list":
            stype = "semi-colon separated list of comma separated lists"
         else:
            stype = str(prop.getType())

         stream.write("### %s (%s)\n" % (prop, stype))
         for line in prop.getDescription().splitlines(False):
            text = "# %s" % (line.expandtabs(3))
            stream.write("%s\n" % text)

         if prop.isExperimental():
            stream.write("### This property is marked as EXPERIMENTAL !!!\n")


      if comments is not None:
         for cline in comments.splitlines():
            stream.write("# %s\n" % cline)
      stream.write("# %s\n" % datetime.datetime.now())

      keys = [key for key in self]

      for key in info:
         if key in keys:
            stream.write("# %s=%s\n" % (key, self[key]))
            keys.remove(key)
      stream.write("\n")

      if properties is not None:
         for propty in properties:
            key = str(propty)
            if propty.isHidden() and not (self[key] is propty.getDefault()):
               continue
            format_comment(propty)

            if key in keys:
               stream.write("%s=%s\n" % (escape(key, True), escape(self[key])))
               keys.remove(key)
            else:
               dflt = propty.getDefault()
               if dflt is None:
                  dflt = ""
               else:
                  dflt = str(dflt)
               stream.write("### Default:\n")
               stream.write("# %s=%s\n" % (escape(key, True), escape(dflt)))
            stream.write("\n")

      stream.write("### Further settings:\n")
      keys.sort()
      for key in keys:
         stream.write("%s=%s\n" % (escape(key, True), escape(self[key])))

   def load(self, stream):
      dollar_variables = []
      buf = ""
      for line in stream:
         buf += line.strip()
         if buf.startswith("!") or buf.startswith("#"):
            buf = ""
            continue
         elif buf == "":
            continue
         elif buf.endswith("\\"):
            bksl = ""
            while buf.endswith("\\\\"):
               bksl += "\\\\"
               buf = buf.rstrip("\\\\")
            if buf.endswith("\\"):
               buf = buf.rstrip("\\")
               buf += bksl
               continue
            else:
               buf += bksl
         # consider the line as terminated now

         # don't do replacements just look for '\\' or
         # one of '=', ':', ' ', '\t', '\f'
         preceeding_backslash = False
         separator_index = -1
         in_brackets = False
         for i in range(len(buf)):
            if (not preceeding_backslash and not in_brackets) and \
                  (buf[i] in ['=', ':', ' ', '\f', '\t']):
               separator_index = i
               break
            if buf[i] == "[":
               in_brackets = True
            elif buf[i] == "]":
               in_brackets = False
            elif buf[i] == '\\':
               preceeding_backslash = not preceeding_backslash
            else:
               preceeding_backslash = False

         if in_brackets:
            raise GolemConfigError("Invalid range specification in '%s'. Missing ']'?" % buf)

         if separator_index < 0:
            key = buf.strip()
            value = ""
         else:
            key = buf[0:separator_index].strip()
            value = buf[separator_index + 1:].strip()

         for dkey, dvalue in reversed(dollar_variables):
            for fmt in ["${%s}", "$(%s)", "$%s"]:
               key = key.replace(fmt % dkey, dvalue)
               value = value.replace(fmt % dkey, dvalue)

         if key.startswith("$"):
            dollar_variables.append( (key[1:], value) )
         else:
            self.setProperty(unescape(key), unescape(value))
         buf = ""

   def __getitem__(self, key):
      return self._getProperty(key)

   def __setitem__(self, key, value):
      self.setProperty(key, value)
   
   def __iter__(self):
      return self.propertyNames()

   def addAll(self, other):
      for name in other:
         self[name] = other[name]
      return self

   def __iadd__(self, other):
      return self.addAll(other)

   def __str__(self):
      res = ""
      for key in self:
         res += "%s=%s\n" % (escape(key, True), escape(self[key]))
      return res

   def _del(self, name):
      del self._map[name]
      # keep plussed and unplussed entries consistent
      if name.startswith("+"):
         del self._map[name[1:]]
      else:
         try:
            del self._map["+" + name]
         except KeyError:
            pass

   def activate_subconfig(self, no):
      changed={}
      for key in self:
        if "[" not in key:
           continue
        pos=key.index("[")
        if no in extractRange(key[pos+1:-1]):
           if key[:pos] in changed.keys() and changed[key[:pos]]!=self[key]:
              raise GolemConfigError("multiple values for option '%s' in subprocess %s: '%s' or '%s'?" \
                    % (key[:pos],no, changed[key[:pos]],self[key]))
           self[key[:pos]]=self[key]
           changed[key[:pos]]=self[key]

def unescape(s):
   buf = [s[i] for i in range(len(s))]
   if "\\" in buf:
      idx_bk = buf.index("\\")
   else:
      idx_bk = -1
   while idx_bk >= 0:
      if idx_bk + 1 < len(buf):
         ch = buf[idx_bk + 1]
         if ch == 'n':
            buf[idx_bk: idx_bk + 2] = "\n"
         elif ch == 'r':
            buf[idx_bk: idx_bk + 2] = "\r"
         elif ch == 'f':
            buf[idx_bk: idx_bk + 2] = "\f"
         elif ch == 't':
            buf[idx_bk: idx_bk + 2] = "\t"
         else:
            del buf[idx_bk: idx_bk + 1]
      if "\\" in buf[idx_bk + 1:]:
         idx_bk = buf.index("\\", idx_bk + 1)
      else:
         idx_bk = -1
   return "".join(buf)

def escape(s, isKey=False):
   escapes = {"\n": "\\n", "\r": "\\r", "\f": "\\f", "\t": "\\t"}
   keyescapes = {"=": "\\=", ":": "\\:", " ": "\\ "}
   buf = s.replace("\\", "\\\\")
   for ch, esc in escapes.items():
      buf = buf.replace(ch, esc)
   if isKey:
      for ch, esc in keyescapes.items():
         buf = buf.replace(ch, esc)
      if buf[0] in ['#', '!']:
         buf = "\\" + buf
   return buf

REQUIRED = 1
OPTIONAL = 0

class ConfigurationException(Exception):
   pass

class Component:
   def __init__(self):
      pass

   def examine(self, hints, **opts):
      pass

   def isInstalled(self):
      return False

   def getInstallationPath(self):
      return []

   def getInstance(self):
      pass

   def undohome(self, name):
      if not name:
         return name
      home = os.getenv("HOME")
      if home is not None:
         if name == home:
            return "${HOME}"
         if name.startswith(home):
            lh = len(home)
            return "${HOME}" + name[lh:]
      return name

   def dohome(self, name):
      if not name:
         return name
      home = os.getenv("HOME")
      return name.replace("${HOME}",home)


   def store(self, conf):
      pass

   def generatePossibleDirs(self, suffix, **opts):
      user_home = gpath.get_homedir()
      golem_dir = gpath.golem_path()
      curdir = os.path.abspath(os.path.curdir)
      result = []

      if suffix == "lib" and "libdir" in opts:
         libdir = opts["libdir"]
         result.append(libdir)

      if suffix == "bin" and "bindir" in opts:
         bindir = opts["bindir"]
         result.append(bindir)

      if "prefix" in opts:
         prefix = opts["prefix"]
         d = os.path.join(prefix, suffix)
         if d not in result:
            result.append(d)

      for path in [ os.path.join(golem_dir, suffix) ]:
         if path not in result:
            result.append(path)

      if suffix == "bin":
         ldp = os.getenv("PATH")
         if ldp:
            for path in ldp.split(os.path.pathsep):
               if path not in result:
                  result.append(path)


      if suffix == "lib":
         ldp = os.getenv("LD_LIBRARY_PATH")
         if ldp:
            for path in ldp.split(os.path.pathsep):
               if path not in result:
                  result.append(path)

      if suffix.startswith("include"):
         # guess include path from LD_LIBRARY_PATH
         ldp = os.getenv("LD_LIBRARY_PATH")
         if ldp:
            for path in ldp.split(os.path.pathsep):
               guess_path= os.path.abspath(os.path.join(path,os.path.pardir,suffix))
               if guess_path not in result:
                  result.append(guess_path)

      for path in [
            os.path.join(golem_dir, suffix),
            os.path.join("/usr", suffix),
            os.path.join("/usr", "local", suffix),
            os.path.join(user_home, suffix),
            os.path.join(user_home, "local", suffix),
            os.path.join(curdir, suffix),
            os.path.join(curdir, "local", suffix)]:
         if path not in result:
            result.append(path)

      return result

class Library(Component):
   def __init__(self, cname, *names):
      self.locations = []
      self.libnames = list(names)
      self.extensions = [".so", ".a", ".dll"]
      self.name = cname

   def examine(self, hints, **opts):
      user_home = gpath.get_homedir()
      golem_dir = gpath.golem_path()
      curdir = os.path.abspath(os.path.curdir)
      dirs = self.generatePossibleDirs("lib")

      if self.name in hints:
         dirs = [hints[self.name]] + dirs

      for dir in dirs:
         found = False
         for libname in self.libnames:
            for ext in self.extensions:
               fname = os.path.join(dir, libname + ext)
               if os.path.exists(fname):
                  self.locations.append(dir)
                  found = True
                  break
            if found:
               break

   def findIncludeDir(self, subpath, incname, hints, *extensions):
      user_home = gpath.get_homedir()
      golem_dir = gpath.golem_path()
      curdir = os.path.abspath(os.path.curdir)
      dirs = self.generatePossibleDirs(os.path.join("include", subpath))

      locations = []

      if self.name in hints:
         dirs = [hints[self.name]] + dirs

      for dir in dirs:
         found = False
         for ext in extensions:
            fname = os.path.join(dir, incname + ext)
            if os.path.exists(fname):
               locations.append(dir)
               found = True
               break
            if found:
               break
      return locations

   def getInstance(self):
      dirs = self.getInstallationPath()
      for libname in self.libnames:
         for ext in self.extensions:
            fname = os.path.join(dirs[0], libname + ext)
            if os.path.exists(fname):
               return fname

   def isInstalled(self):
      return len(self.locations) > 0

   def getInstallationPath(self):
      return self.locations[:]

class Program(Component):
   def __init__(self, cname, *names):
      self.locations = []
      self.prognames = list(names)
      self.extensions = ["", ".exe", ".com"]
      self.name = cname

   def examine(self, hints):
      user_home = gpath.get_homedir()
      golem_dir = gpath.golem_path()
      curdir = os.path.abspath(os.path.curdir)
      dirs = self.generatePossibleDirs("bin")

      if self.name in hints:
         hint = hints[self.name]
         if os.path.exists(hint):
            self.locations.append(hint)

      for progname in self.prognames:
         for ext in self.extensions:
            for dir in dirs:
               fname = os.path.join(dir, progname + ext)
               if os.path.exists(fname):
                  if dir not in self.locations:
                     self.locations.append(fname)

   def getInstance(self,conf=None):
      if len(self.locations) > 0:
         return self.locations[0]

   def isInstalled(self):
      return len(self.locations) > 0

   def getInstallationPath(self):
      return self.locations[:]

class LoopTools(Library):
   def __init__(self):
      Library.__init__(self, "LoopTools", "libooptools")

   def store(self, conf):
      paths = self.getInstallationPath()
      if len(paths) == 0:
         return

      path = self.undohome(paths[0])

      if "+zzz.extensions" in conf:
         conf["+zzz.extensions"] += ", looptools"
      else:
         conf["+zzz.extensions"] = "looptools"

      conf["+looptools.ldflags"] = "-L%s -looptools" % path

class AVH_OneLoop(Library):
   def __init__(self):
      Library.__init__(self, "avholo", "libavh_olo")

   def store(self, conf):
      paths = self.getInstallationPath()
      if len(paths) == 0:
         return

      path = self.undohome(paths[0])

      if "+zzz.extensions" in conf:
         conf["+zzz.extensions"] += ", avh_olo"
      else:
         conf["+zzz.extensions"] = "avh_olo"

      conf["+avh_olo.ldflags"] = "-L%s -lavh_olo" % path

class QCDLoop(Library):
   def __init__(self):
      Library.__init__(self, "QCDLoop", "libqcdloop")

   def store(self, conf):
      paths = self.getInstallationPath()
      if len(paths) == 0:
         return

      path = self.undohome(paths[0])

      if "+zzz.extensions" in conf:
         conf["+zzz.extensions"] += ", qcdloop"
      else:
         conf["+zzz.extensions"] = "qcdloop"

      conf["+qcdloop.ldflags"] = "-L%s -lqcdloop" % path

class Ninja(Library):
   def __init__(self):
      Library.__init__(self, "Ninja", "libninja")

   def examine(self, hints):
      Library.examine(self, hints)
      if len(self.locations) > 0:
         self.incdirs = self.findIncludeDir("gosam-contrib", "ninjago_module", hints,
               ".mod") or self.findIncludeDir("ninja-contrib", "ninjago_module", hints,
               ".mod")
         if len(self.incdirs) == 0:
            self.locations = []

   def store(self, conf):
      paths = self.getInstallationPath()
      if len(paths) == 0:
         return

      path = self.undohome(paths[0])
      incd = self.undohome(self.incdirs[0])

      if "+installed.extensions" in conf:
         conf["+installed.extensions"] += ", ninja"
      else:
         conf["+installed.extensions"] = "ninja"

      conf["ninja.fcflags"] = "-I%s" % incd
      conf["ninja.ldflags"] = "-L%s -lninja" % path



class Samurai(Library):
   def __init__(self):
      Library.__init__(self, "Samurai", "libsamurai")

   def examine(self, hints):
      Library.examine(self, hints)
      if len(self.locations) > 0:
         self.incdirs = self.findIncludeDir("gosam-contrib", "msamurai", hints,
               ".mod") or self.findIncludeDir("samurai", "msamurai", hints,
               ".mod")
         if len(self.incdirs) == 0:
            self.locations = []

   def store(self, conf):
      paths = self.getInstallationPath()
      if len(paths) == 0:
         return

      path = self.undohome(paths[0])
      incd = self.undohome(self.incdirs[0])

      if "+installed.extensions" in conf:
         conf["+installed.extensions"] += ", samurai"
      else:
         conf["+installed.extensions"] = "samurai"

      conf["samurai.fcflags"] = "-I%s" % incd
      conf["samurai.ldflags"] = "-L%s -lsamurai" % path

class PJFry(Library):
   def __init__(self):
      Library.__init__(self, "PJFry", "libpjfry")

   def store(self, conf):
      paths = self.getInstallationPath()
      if len(paths) == 0:
         return

      path = self.undohome(paths[0])

      if "+installed.extensions" in conf:
         conf["+installed.extensions"] += ", pjfry"
      else:
         conf["+installed.extensions"] = "pjfry"

      conf["+pjfry.ldflags"] = "-L%s -lpjfry" % path

class Golem95(Library):
   def __init__(self):
      Library.__init__(self, "Golem", "libgolem")

   def examine(self, hints):
      Library.examine(self, hints)
      if len(self.locations) > 0:
         self.incdirs = self.findIncludeDir("gosam-contrib", "parametre", hints,
               ".mod") or self.findIncludeDir("golem95", "parametre", hints,
               ".mod")
         if len(self.incdirs) == 0:
            self.locations = []

   def store(self, conf):
      paths = self.getInstallationPath()
      if len(paths) == 0:
         return

      path = self.undohome(paths[0])
      incd = self.undohome(self.incdirs[0])

      if "+installed.extensions" in conf:
         conf["+installed.extensions"] += ", golem95"
      else:
         conf["+installed.extensions"] = "golem95"

      conf["golem95.fcflags"] = "-I%s" % incd
      conf["golem95.ldflags"] = "-L%s -lgolem" % path

class Form(Program):
   def __init__(self):
      Program.__init__(self, "Form",
            "tformi", "formi", "tvormi", "vormi",
            "tform", "form", "tvorm", "vorm",
            "tform3", "form3", "tvorm3", "vorm3")

   def examine(self, hints):
      Program.examine(self, hints)
      executable = self.getInstance()

      if executable is None:
         return

      try:
         pipe = subprocess.Popen(executable,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         (stdout, stderr) = pipe.communicate()

         stdout_lines = stdout.splitlines()

         firstline=True
         for line in stdout_lines:
            lline = line.lower().strip()
            if "version" in lline:
               i = lline.index("version")
               j = lline.index("(")
               lline = lline[i+len("version"):j]
               self.version = [int(re.sub("[^0-9]", "", s))
                     for s in lline.split(".")]
            elif firstline and "form" in lline:
               # form version 4.0 prints
               # "FORM 4.0 (Jun 11 2012) 64-bits"
               try:
                  lline = lline.split(" ")[1]
               except IndexError:
                  pass
               self.version = [int(re.sub("[^0-9]", "", s))
                     for s in lline.split(".")]
            firstline=False

      except OSError:
         raise ConfigurationException(
               "Could not run FORM (%s) properly" % executable)

   def checkCompatibility(self,p):
      #if not testFormCompatibility(p):
      #   return False
      return True

   def store(self, conf):
      conf["form.bin"] = self.undohome(self.getInstance(check=True))
      #if version_compare(self.version, [4,0]) >= 0:
      #   conf["+form.extensions"] = "topolynomial"

   def getInstance(self,check=False):
      if check:
         for p in self.locations:
            print "# ~~~ " + p + " usable with GoSam? ... ",
            if self.checkCompatibility(p):
               print "Yes"
               return p
            else:
               print "No"
         print("==> Configuration failed:")
         print("    FORM version is not compatible with GoSam.")
         sys.exit(1)
      elif len(self.locations) > 0:
         return self.locations[0]

class Yaml(Component):
   def __init__(self):
      Component.__init__(self)
      self.installed = False

   def examine(self, hints, **opts):
      try:
         import yaml
         self.installed = True
      except ImportError:
         return

   def isInstalled(self):
      return self.installed

   def getInstallationPath(self):
      try:
         import yaml
         return [str(yaml.__file__)]
      except ImportError:
         return []

class Reduze(Program):
   def __init__(self):
      Program.__init__(self, "Reduze", "reduze")

   def examine(self, hints):
      Program.examine(self, hints)
      executable = self.getInstance()

      if executable is None:
         return

      try:
         pipe = subprocess.Popen(executable,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
         (stdout, stderr) = pipe.communicate()

      except OSError:
         raise ConfigurationException(
              "Could not run Reduze (%s) properly" % executable)

   def checkCompatibility(self,p):
      pipe = subprocess.Popen([p, "-v"],
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (stdout, stderr) = pipe.communicate()

      stdout_lines = stdout.splitlines()

      # first line
      lline = stdout_lines[0].lower().strip()
      if "reduze" not in lline:
         # reduze version 2.1 prints
         # "Reduze 2.1.x-297-g01c1202 (non-MPI build)"
         return False
      if not testReduzeCompatibility(p):
         return False
      return True

   def getInstance(self,check=False):
      if check:
         for p in self.locations:
            print "# ~~~ " + p + " usable with GoSam? ... ",
            if self.checkCompatibility(p):
               print "Yes"
               return p
            else:
               print "No"
               return
      elif len(self.locations) > 0:
         return self.locations[0]

   def store(self, conf):
      conf["reduze.bin"] = self.undohome(self.getInstance(check=True))

class SymPy(Component):
   def __init__(self):
      Component.__init__(self)
      self.installed = False

   def examine(self, hints, **opts):
      try:
         import sympy
         self.installed = True
      except ImportError:
         return

   def isInstalled(self):
      return self.installed

   def getInstallationPath(self):
      try:
         import sympy
         return [str(sympy.__file__)]
      except ImportError:
         return []

class Fortran(Program):
   def __init__(self):
      Program.__init__(self, "Fortran",
            "ifort", "f95i",
            "f95", "f95n",
            "lfc", "lf95", "f95f",
            "xlf95", "xlf90", "xlf",
            "gfortran", "g95",
            "f95", "f90", "frt", "pgf", "pgf90", "pghpf",
            "fort90", "fl64", "fl32",
            "pgf77", "g77", "fort77", "f77", "af77", "f2c")

   def examine(self, hints):
      fc = os.getenv("FC")
      files=["gosam.conf"]
      directories = [gpath.golem_path(), gpath.gosam_contrib_path()]
      for dir in directories:
         for file in files:
            full_name = os.path.join(dir, file)
            if os.path.exists(full_name):
               for l in open(full_name):
                  m=re.match("^\s*fc.bin=\s*([^#]+)\s*(#.*)?$",l.strip())
                  if m:
                     path, prog = os.path.split(m.group(1))
                     path=path.strip()
                     prog=prog.strip()
                     if len(prog) > 0:
                        self.prognames = [prog] + self.prognames
                     if len(path) > 0 and self.name not in hints:
                        hints = hints.copy()
                        hints[self.name] = path
                     if len(prog)  > 0 and len(path) > 0 and os.path.exists(m.group(1)):
                        self.locations.append(m.group(1))

      if fc is not None:
         path, prog = os.path.split(fc)
         if len(prog) > 0:
            self.prognames = [prog] + self.prognames
         if len(path) > 0 and self.name not in hints:
            hints = hints.copy()
            hints[self.name] = path

      Program.examine(self, hints)

   def checkCompatibility(self,p,conf):
      p=self.dohome(p)
      if "golem95.fcflags" in conf:
         if not testCompilerLibCompatibility(p,'parametre',self.dohome(conf["golem95.fcflags"])):
            return False
      if "samurai.fcflags" in conf:
         if not testCompilerLibCompatibility(p,'msamurai',self.dohome(conf["samurai.fcflags"])):
            return False
      return True

   def getInstance(self,conf=None):
      if conf:
         for p in self.locations:
            print "# ~~~ " + p + " usable with installed Golem95/Samurai? ... ",
            if self.checkCompatibility(p,conf):
               print "Yes"
               return p
            else:
              print "No"
         print("==> Configuration failed:")
         print("    Libraries not with examined compiler created.")
         sys.exit(1)
      elif len(self.locations) > 0:
         return self.locations[0]


   def store(self, conf):
      conf["fc.bin"] = self.undohome(self.getInstance(conf=conf))

class QGraf(Program):
   def __init__(self):
      Program.__init__(self, "QGraf", "qgraf")

   def store(self, conf):
      conf["qgraf.bin"] = self.undohome(self.getInstance())

class Java(Program):
   def __init__(self):
      Program.__init__(self, "Java", "java")

   def store(self, conf):
      haggies_jar = self.undohome(gpath.golem_path("haggies", "haggies.jar"))
      java = self.undohome(self.getInstance())
      conf["haggies.bin"] = "%s -jar %s" % \
            (java, haggies_jar)

class Configurator:
   def __init__(self, hints, **components):

      self.installed_components = []

      not_found = []

      items=components.items()

      def preferLibKey(x):
         if x[0]=="Golem95":
            return "  " + x[0]
         if x[0]=="Samurai":
            return " " + x[0]
         return x[0]

      # search first for Golem95/Samurai to find a compatible compiler
      items.sort(key=preferLibKey)

      for name, required in items:
         if name not in globals():
            raise ConfigurationException("Name '%s' not known" % name)

         cls = globals()[name]

         if type(cls) != type(Component):
            raise ConfigurationException(
                  "Name '%s' does not denote a class" % name)
         if not issubclass(cls, Component):
            raise ConfigurationException(
                  "Class '%s' is not a subclass of 'Component'" % name)
                  
         component = cls()

         self.message("Searching for %s ..." % name)
         component.examine(hints)

         if not component.isInstalled():
            if required == REQUIRED:
               not_found.append(name)
               self.message("Required component %s is not installed." % name)
            else:
               self.message("Component %s is not installed." % name)
         else:
            paths = component.getInstallationPath()
            l = len(paths)
            if l == 1:
               self.message("Component %s has been found." % name)
               self.message("     %s" % (paths[0]))
            else:
               self.message("Component %s has been found in %d places." %
                     (name, l))
               for i, path in enumerate(paths):
                  self.message("#%2d: %s" % (i+1, path))
            self.installed_components.append(component)

      if len(not_found) > 0:
         self.fail("Required component(s): %s not installed." % ', '.join(not_found))

   def message(self, message):
      print("# ~~~ " + message)

   def fail(self, message):
      print("==> Configuration failed:")
      print("    " + message)
      sys.exit(1)

   def store(self, conf):
      for component in self.installed_components:
         component.store(conf)

def version_compare(v1, v2):
   l1 = len(v1)
   l2 = len(v2)
   if l1 < l2:
      v1c = v1 + [0] * (l2-l1)
      v2c = v2[:]
   else:
      v1c = v1[:]
      v2c = v2 + [0] * (l1-l2)

   for p1, p2 in zip(v1c, v2c):
      if p1 > p2:
         return 1
      elif p1 < p2:
         return -1
   return 0

def levenshtein(str1, str2, case_sensitive=False):
	l1 = len(str1)
	l2 = len(str2)

	if case_sensitive:
		if l1 < l2:
			if l1 == 0:
				return l2
			s2 = str1
			s1 = str2
		else:
			if l2 == 0:
				return l1
			s1 = str1
			s2 = str2
	else:
		if l1 < l2:
			if l1 == 0:
				return l2
			s2 = str1.lower()
			s1 = str2.lower()
		else:
			if l2 == 0:
				return l1
			s1 = str1.lower()
			s2 = str2.lower()


	previous_row = xrange(len(s2) + 1)
	for i, c1 in enumerate(s1):
		current_row = [i + 1]
		for j, c2 in enumerate(s2):
			# j+1 instead of j since previous_row and current_row
			# are one character longer
			insertions = previous_row[j + 1] + 1
			# than s2
			deletions = current_row[j] + 1
			substitutions = previous_row[j] + (c1 != c2)
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row

	return previous_row[-1]

def testCompilerLibCompatibility (compiler,lib,flags):
   cur_path=os.getcwd()
   try:
      tmp_dir=tempfile.mkdtemp()
      os.chdir(tmp_dir)
      with open("program.f90","w") as f:
         f.write("program test\n")
         f.write("use %s\n" % lib)
         f.write("end program test\n")
      p=subprocess.Popen([compiler,"program.f90","-c"] + shlex.split(flags),
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      (stdout, stderr) = p.communicate()
      return p.returncode==0
   finally:
      os.chdir(cur_path)
      try:
         shutil.rmtree(tmp_dir) # delete directory
      except OSError, e:
         if e.errno != 2: # no such file or directory
            raise

def testFormCompatibility(executable):
   cur_path=os.getcwd()
   try:
      tmp_dir = tempfile.mkdtemp()
      os.chdir(tmp_dir)

      # Check that required newer FORM features work

      #  1) PolyRatFun rat(divergence,x)
      with open("divergence.frm","w") as f:
         f.write("CF rat;\n")
         f.write("S x;\n")
         f.write("PolyRatFun rat(divergence, x);\n")
         f.write("L e = rat(x+1,x^2);\n")
         f.write(".sort\n")
         f.write("#Write <divergence.out> \"%e\",e\n")
         f.write(".end\n")
      p = subprocess.Popen([executable, "divergence.frm"],
                           stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (stdout, stderr) = p.communicate()
      divergence_output_expected =  "rat(1,x^2);"
      with open("divergence.out","r") as f:
         divergence_output = f.read()
      divergence_output = divergence_output.strip()
      if divergence_output != divergence_output_expected:
         return False

      return p.returncode == 0

   except:
      return False

   finally:
      os.chdir(cur_path)
      try:
         shutil.rmtree(tmp_dir)  # delete directory
      except OSError, e:
         if e.errno != 2:  # no such file or directory
            raise

def testReduzeCompatibility(executable):
   try:
      from yaml import load, dump
   except ImportError:
      return False

   cur_path=os.getcwd()
   try:
      # setup a very basic reduction and match a basic diagram, then check all output is as we expect
      tmp_dir=tempfile.mkdtemp()
      os.chdir(tmp_dir)
      os.mkdir(tmp_dir+'/config')
      os.chdir(tmp_dir+'/config')

      with open("integralfamilies.yaml","w") as f:
         f.write("integralfamilies:\n")
         f.write("  - name: triangle\n")
         f.write("    loop_momenta: [k1]\n")
         f.write("    propagators:\n")
         f.write("      - [\"k1\",0]\n")
         f.write("      - [\"k1+p1\",0]\n")
         f.write("      - [\"k1+p1+p2\",0]\n")

      with open("kinematics.yaml","w") as f:
         f.write("kinematics:\n")
         f.write("  incoming_momenta: [p1, p2, p3]\n")
         f.write("  outgoing_momenta: []\n")
         f.write("  momentum_conservation: [p3, -p1 -p2]\n")
         f.write("  kinematic_invariants:\n")
         f.write("    - [s,  2]\n")
         f.write("    - [m,  1]\n")
         f.write("  scalarproduct_rules:\n")
         f.write("    - [[p1,p1],\"0\"]\n")
         f.write("    - [[p2,p2],\"0\"]\n")
         f.write("    - [[p1+p2,p1+p2],\"m^2\"]\n")
         f.write("  symbol_to_replace_by_one: s\n")

      os.chdir(tmp_dir)

      with open("diagrams.yaml","w") as f:
         f.write("---\n")
         f.write("qgraf_globals:\n")
         f.write("  - \"in=phi[p1],phi[p2],phim[p3];\"\n")
         f.write("  - \"out=;\"\n")
         f.write("  - \"loops=1;\"\n")
         f.write("  - \"loop_momentum = k;\"\n")
         f.write("\n")
         f.write("---\n")
         f.write("diagram:\n")
         f.write("  name: 5\n") # Note: diagram is not diagram 1, this checks if Reduze is renumbering diagrams (it should not)
         f.write("  external_legs:\n")
         f.write("    - [ [-1, 1], in-phi: [[phi, 1, +1, p1, 0, -1]] ]\n")
         f.write("    - [ [-2, 2], in-phi: [[phi, 1, +1, p2, 0, -3]] ]\n")
         f.write("    - [ [-3, 3], in-phim: [[phim, 1, +1, p3, m, -5]] ]\n")
         f.write("  propagators:\n")
         f.write("    - [ [2, 1], phi_phi: [[phi, 3, +1, -k1, 0, 1], [phi, 3, +1, k1, 0, 2]] ]\n")
         f.write("    - [ [3, 1], phi_phi: [[phi, 3, +1, k1-p1, 0, 3], [phi, 3, +1, -k1+p1, 0, 4]] ]\n")
         f.write("    - [ [3, 2], phi_phi: [[phi, 3, +1, -k1-p2, 0, 5], [phi, 3, +1, k1+p2, 0, 6]] ]\n")
         f.write("  vertices:\n")
         f.write(
            "    - [ phi_phi_phi: [[phi, 1, +1, p1, 0, -1], [phi, 3, +1, -k1, 0, 1], [phi, 3, +1, k1-p1, 0, 3]] ]\n")
         f.write(
            "    - [ phi_phi_phi: [[phi, 1, +1, p2, 0, -3], [phi, 3, +1, k1, 0, 2], [phi, 3, +1, -k1-p2, 0, 5]] ]\n")
         f.write(
            "    - [ phi_phi_phim: [[phi, 3, +1, -k1+p1, 0, 4], [phi, 3, +1, k1+p2, 0, 6], [phim, 1, +1, p3, m, -5]] ]\n")
         f.write("  symmetry_factor: 1\n")
         f.write("  num_legs_in: 3\n")
         f.write("  num_legs_out: 0\n")
         f.write("  num_loops: 1\n")
         f.write("  num_propagators: 3\n")
         f.write("  num_vertices: 3\n")

      with open("job.yaml","w") as f:
         f.write("jobs:\n")
         f.write("  - setup_sector_mappings:\n")
         f.write("      construct_minimal_graphs: true\n")
         f.write("  - setup_sector_mappings_alt:\n")
         f.write("      source_sectors:\n")
         f.write("        select_all: true\n")
         f.write("      find_sector_symmetries: true\n")
         f.write("  - find_diagram_shifts:\n")
         f.write("      qgraf_file: \"diagrams.yaml\"\n")
         f.write("      output_file: \"diagrams.match.yaml\"\n")
         f.write("      info_file_form: \"diagrams.match.inc\"\n")
      p=subprocess.Popen([executable,"job.yaml"],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      (stdout, stderr) = p.communicate()

      # Check if the diagram was matched and if the name of the diagram was correctly maintained
      with open("diagrams.match.yaml","r") as f:
         diagrams_match = load(f)
      if diagrams_match['diagram']['name'] != 5:
         return False

      # Check sectormappings folder was created by Reduze
      os.chdir(tmp_dir+"/sectormappings")

      # Check crossings.yaml has name key (created by Reduze)
      # Note: can not check entire crossings.yaml as rules_momenta are not written deterministically
      crossings_output_expected_name = "x12"
      with open("crossings.yaml","r") as f:
         crossings_output = load(f)
      if crossings_output['crossings'][0]['ordered_crossings'][0]['name'] != crossings_output_expected_name:
         return False

      # Check crss.yaml is correct (created by Reduze)
      crss_output_expected = \
         [[[
            {'sector_selection':
                {'deselect': [],
                 'select_all': False,
                 'deselect_recursively': [],
                 'deselect_independents': [],
                 'deselect_graphless': False,
                 't_restriction': [-1, -1],
                 'select_recursively': [['triangle', 7]],
                 'select': []
                 }
             }
         ]]]
      with open("crss.yaml","r") as f:
         crss_output = load(f)
      if crss_output != crss_output_expected:
         return False

      return p.returncode==0

   except:
      return False

   finally:
      os.chdir(cur_path)
      try:
         shutil.rmtree(tmp_dir) # delete directory
      except OSError, e:
         if e.errno != 2: # no such file or directory
            raise

def extractRange(s,minval=0,maxval=999):
   """ extract ranges from a string and returns a set
       Examples:

      >>> sorted(extractRange("2-3,5-8"))
      [2, 3, 5, 6, 7, 8]
      >>> sorted(extractRange("2-8,!4,!10"))
      [2, 3, 5, 6, 7, 8]
      >>> sorted(extractRange("2-8, !4-6"))
      [2, 3, 7, 8]
      >>> sorted(extractRange("1-10", maxval=5))
      [1, 2, 3, 4, 5]
      >>> sorted(extractRange("2-", 3, 5))
      [3, 4, 5]
      >>> sorted(extractRange("-2"))
      [0, 1, 2]
      >>> extractRange("!4")
      set([])
      >>> sorted(extractRange("-")) == range(0,1000)
      True

   """
   s=s.replace(" ",",")
   ranges = [x.split("-") for x in s.split(",")]
   res=set()
   for r in ranges:
      if r==['']:
         continue
      elif len(r)==1:
         if r[0][0] in "^!":
            res.discard(int(r[0][1:]))
         else:
            res.add(int(r[0]))
      elif len(r)==2:
         remove = r[0] and (r[0][0] in "^!")
         if remove:
            r[0]=r[0][1:]
         start = minval if r[0]=='' else int(r[0])
         end = maxval if r[1]=='' else int(r[1])
         if remove:
            res=set(filter(lambda x: x<start or x>end ,res))
         else:
            res.update(range(start,end+1))
      else: # len(r)>2
         raise ValueError("Invalid range: %s in '%s'" % (r,s))
   res=set(filter(lambda x: x>=minval and x<=maxval,res))
   return res

def split_qgrafPower(power):
   """
   >>> split_qgrafPower('QCD,2,0,QED,3,3')
   [['QCD', 2, 0], ['QED', 3, 3]]
   >>> split_qgrafPower('QCD,2,QED,3')
   [['QCD', 2], ['QED', 3]]
   >>> split_qgrafPower('QCD,2,3,QED,3')
   [['QCD', 2, 3], ['QED', 3, 3]]
   >>> split_qgrafPower('QCD,2')
   [['QCD', 2]]
   >>> split_qgrafPower('QED,3,4')
   [['QED', 3, 4]]
   >>> split_qgrafPower('QCD,2,3,4,QED,3,NP,1')
   [['QCD', 2, 3, 4], ['QED', 3, 3, 3], ['NP', 1, 1, 1]]
   >>> split_qgrafPower('QED,3,4,QED,3,4')
   Traceback (most recent call last):
    ...
   ConfigurationException: Coupling 'QED' repeated
   """
   if type(power)==list:
      return power
   assert(type(power)==str)
   min_length=0
   orders=[]
   couplings=set()
   l=re.split(',|;',power)
   current_coupling=[]
   for i in l + [""]:
      if str(i).isdigit() or str(i).lower()=="none":
         assert(current_coupling)
         current_coupling.append(i)
      else:
         current_len=len(current_coupling)-1
         if current_len<0:
            current_len=0
         if current_len<min_length and current_coupling:
            if current_len>0:
                current_coupling.extend([current_coupling[-1]]*(min_length-current_len))
            else:
                current_coupling.extend([0]*(min_length-current_len))
         if current_len and current_coupling:
             orders.append(current_coupling)
         if min_length==0:
            min_length=current_len
         if i:
             current_coupling=[i]
             if i in couplings:
                 raise ConfigurationException("Coupling '%s' repeated" % i)
             couplings.add(i)
         else:
             current_coupling=[]
   return orders
