# vim: ts=3:sw=3:expandtab
from golem.model.particle import Particle
import golem.util.tools
from math import sqrt

model_name = "Scalar QED in Feynman Gauge (photon, charged leptons, phi+-; EW counterterm insertions)"
#---#[ particles:
# The last argument of Particle is the electric charge of the field in units
# of e/3, with its own sign (e+ = +3, phi- = -3): see charge_convention.
particles = {
   'A':       Particle('A',        2,      0,  1, None,    '0',       22,     0.0),

   'em':      Particle('em',       1,   'me',  1, 'ep',    '0',       11,    -3.0),
   'ep':      Particle('ep',      -1,   'me',  1, 'em',    '0',      -11,     3.0),

   'mum':     Particle('mum',      1,  'mmu',  1, 'mup',   '0',       13,    -3.0),
   'mup':     Particle('mup',     -1,  'mmu',  1, 'mum',   '0',      -13,     3.0),

   'taum':    Particle('taum',     1, 'mtau',  1, 'taup', 'wtau',     15,    -3.0),
   'taup':    Particle('taup',    -1, 'mtau',  1, 'taum', 'wtau',    -15,     3.0),

   'phim':    Particle('phim',     0,   'mW',  1, 'phip',    '0',   -251,    -3.0),
   'phip':    Particle('phip',     0,   'mW',  1, 'phim',    '0',    251,     3.0),
   'Cx':      Particle('Cx',      0,      0,  1, None, '0', 0,     0.0)   
}
#---#] particles:
# Charges are the physical charges of each field (see above), not the
# Standard-Model convention of smdiag_complex_ct (e+ = -3, phi+- = +-1, read
# through hard-coded PDG rules); used by the QED infrared operator
# (ir_subtraction) and by the conversion to CDR.
charge_convention = "physical"
#---#[ mnemonics:
mnemonics = {
   'e+': particles['ep'],
   'e-': particles['em'],
   'positron': particles['ep'],
   'electron': particles['em'],
   'mu+': particles['mup'],
   'mu-': particles['mum'],
   'tau+': particles['taup'],
   'tau-': particles['taum'],
   'photon': particles['A'],
   'gamma': particles['A'],
   'phi-': particles['phim'],
   'phi+': particles['phip']
}
#---#] mnemonics:
#---#[ parameters:
# Parameters from the Particle Data Booklet 2008
# see http://pdg.lbl.gov/
# <---- indicates deviations from the PDG table
# Masses in GeV
parameters = {
   'NC': '3.0',
   'gs': '1.0', # <----
   'me': '0.000510998910',
   'mmu': '0.105658367',
   'mtau': '1.77684',
   # Widths <----
   'wZ': '2.4952',
   'wW': '2.124',
   'wtau': '0.0',
   # Number of flavours, not really a model parameter
   # but needed:
   'Nf': '5.0',
   'Nfgen': '-1.0'
}
#---#] parameters:
#---#[ latex_parameters:
latex_parameters = {
   'NC': 'N_C',
   'e': 'e',
   'alpha': "\\alpha",
   'GF': "G_F",
   'gs': 'g_s',
   'me': 'm_e',
   'mmu': 'm_\\mu',
   'mtau': 'm_\\tau',
   'mZ': 'm_Z', 'mW': 'm_W',
   'sw': 's_w',
   'wZ': '\\Gamma_Z', 'wW': '\\Gamma_W',
   'wtau': '\\Gamma_\\tau',
   'Nf': 'N_f',
   'Nfgen': 'N_f^{gen}',
   'NA': 'N_A',
   'Nfrat': 'N_f^{rat}',
   'gAAPP': 'g_{\\gamma\\gamma\\phi^+\\phi^-}',
   'gAPP': 'g_{\\gamma\\phi^+\\phi^-}'
}
#---#] latex_parameters:
#---#[ functions:
functions = {
   'NA': 'NC*NC-1',
   'gAAPP': '2',
   'gAPP': '-1',

   'Nfrat': 'if(Nfgen,Nf/Nfgen,1)'
}
#---#] functions:
#---#[ types:
types = {
   'NC': 'R', 'gs': 'R',
   'me': 'R', 'mmu': 'R', 'mtau': 'R',
   'wtau': 'R',
   'wZ': 'R', 'wW': 'R',
   'Nf': 'R', 'Nfgen': 'R', 'NA': 'R',
   'Nfrat': 'R',
   'gAAPP': 'C',
   'gAPP': 'R'
}
#---#] types:
#---#[ latex_names:
latex_names = {
   'A':       '\\gamma',

   'em':      'e^-',
   'ep':      'e^+',

   'mum':     '\\mu^-',
   'mup':     '\\mu^+',

   'taum':    '\\tau^-',
   'taup':    '\\tau^+',

   'phim':    '\phi^-',
   'phip':    '\phi^+',
   'Cx': 'Cx'
}
#---#] latex_names:
#---#[ line_styles:
line_styles = {
   'A':       'photon',

   'em':      'fermion',
   'ep':      'fermion',

   'mum':     'fermion',
   'mup':     'fermion',

   'taum':    'fermion',
   'taup':    'fermion',

   'phim':    'scalar',
   'phip':    'scalar',
}
#---#] line_styles:
#---#[ slha_locations:
slha_locations = {
   'me': ('MASS', [11]),
   'mmu': ('MASS', [13]),
   'mtau': ('MASS', [15]),
   'mZ': ('MASS', [23]),
   'mW': ('MASS', [24]),
   'wtau': ('DECAY', [15]),
   'wZ': ('DECAY', [23]),
   'wW': ('DECAY', [24]),
}
#---#] slha_locations:

#---#[ def init_ew:
def init_ew(e_one=False,**options):
   """
   Produce entries in parameters and functions starting from the
   given initial values.

   Reference Formulae:

   GF / sqrt(2) = alpha * pi / 2 / mW^2 / sw^2
                = e^2 / 8 / mW^2 / sw^2
        sw^2 = 1 - mW^2 / mZ^2

   We have access to the "user_choice" from options
   and if any parameters were specified we determine
   gosam_choice from the input parameters

   Then we compare...

   """
   global parameters, functions, types

   keys = set(options.keys())
   for key in keys:
      parameters[key] = str(options[key])
      types[key] = "R"
   if keys == set(["GF", "mW", "mZ"]):
      # mW, mZ --> sw
      functions["sw"] = "csqrt(1-(mW*mW-i_*mW*wW)/(mZ*mZ-i_*mZ*wZ))"
      types["sw"] = "C"
      # GF, mW, sw --> e
      functions["e"] = "mW*sw*sqrt(8*GF/sqrt(2))"
      types["e"] = "C"
   elif keys == set(["alpha", "mW", "mZ"]):
      # alpha --> e
      functions["e"] = "sqrt(4*pi*alpha))"
      types["e"] = "R"
      # mW, mZ --> sw
      functions["sw"] = "csqrt(1-(mW*mW-i_*mW*wW)/(mZ*mZ-i_*mZ*wZ))"
      types["sw"] = "C"
   elif keys == set(["alpha", "sw", "mZ"]):
      # alpha --> e
      functions["e"] = "sqrt(4*pi*alpha))"
      types["e"] = "R"
      # sw, mZ --> mW
      functions["mW"] = "mZ*sqrt(1-real(sw*sw))"
      types["mW"] = "R"
   elif keys == set(["alpha", "sw", "GF"]):
      # alpha --> e
      functions["e"] = "sqrt(4*pi*alpha))"
      types["e"] = "R"
      # GF, sw, alpha --> mW
      functions["mW"] = "sqrt(alpha*pi/sqrt(2)/GF) / real(sw)"
      types['mW'] = 'R'
      # mW, sw --> mZ
      functions["mZ"] = "mW / sqrt(1-real(sw*sw))"
      types["mZ"] = "R"
   elif keys == set(["alpha", "GF", "mZ"]):
      # alpha --> e
      functions["e"] = "sqrt(4*pi*alpha))"
      types["e"] = "R"
      # GF, mZ, alpha --> mW
      functions["mW"] = "sqrt(mZ*mZ/2+sqrt(mZ*mZ*mZ*mZ/4-pi*alpha*mZ*mZ/sqrt(2)/GF))"
      types["mW"] = "R"
      # mW, mZ --> sw
      functions["sw"] = "csqrt(1-(mW*mW-i_*mW*wW)/(mZ*mZ-i_*mZ*wZ)"
      types["sw"] = "C"
   elif keys == set(["e", "mW", "mZ"]):
      # mW, mZ --> sw
      functions["sw"] = "csqrt(1-(mW*mW-i_*mW*wW)/(mZ*mZ-i_*mZ*wZ)"
      types["sw"] = "C"
   elif keys == set(["e", "sw", "mZ"]):
      gosam_choice = 7
      # mZ, sw --> mW
      functions["mW"] = "mZ*sqrt(1-real(sw*sw))"
      types["mW"] = "R"
   elif keys == set(["e", "sw", "GF"]):
      # e, sw, GF --> mW
      functions["mW"] = "e/2/sw/sqrt(sqrt(2)*GF)"
      types["mW"] = "R"
      # mW, sw --> mZ
      functions["mZ"] = "mW / sqrt(1-real(sw*sw))"
      types["mZ"] = "R"
   elif keys == set(["e", "sw", "GF", "mZ", "mW", "alpha"]):
      for dummy in ["e", "sw", "GF", "mZ", "mW", "alpha"]:
         #   parameters[dummy] = '0.0'
         functions['%sf' % dummy ] = dummy
         types[dummy] = "R"
         types['%sf' % dummy] = "R"
         types["sw"] = "C"
         types["swf"]= "C"
         #try:
         #   del slha_locations[dummy]
         #except:
         #   continue
   else:
      raise Exception("Invalid EW Scheme.")
   if e_one:
      if "e" in functions.keys():
         del functions["e"]
#---#] def init_ew:

#---#[ def ew_gosam_choice:

def ew_gosam_choice(keys, size_set):
   """
      Returns GoSam's prefered choice of electroweak scheme based on the
      input parameters in "keys"

      The values are 1-8, with the special results
       0 : no parameters specified
      -1 : less than 3 parameters specified (so a unique choice
           cannot be made)
      -2 : more than 3 parameters specified (so a unique choice
           cannot be made)
   """
   # is it one of these distince choices?
   if keys == set(["GF", "mW", "mZ"]):
      gosam_choice = 1
   elif keys == set(["alpha", "mW", "mZ"]):
      gosam_choice = 2
   elif keys == set(["alpha", "sw", "mZ"]):
      gosam_choice = 3
   elif keys == set(["alpha", "sw", "GF"]):
      gosam_choice = 4
   elif keys == set(["alpha", "GF", "mZ"]):
      gosam_choice = 5
   elif keys == set(["e", "mW", "mZ"]):
      gosam_choice = 6
   elif keys == set(["e", "sw", "mZ"]):
      gosam_choice = 7
   elif keys == set(["e", "sw", "GF"]):
      gosam_choice = 8
   elif size_set == 0:
      gosam_choice = 0
   elif size_set < 3 and size_set > 0:
      gosam_choice = -1
   elif size_set > 3:
      gosam_choice = -2
   else:
      raise Exception("Invalid EW Scheme input by user.")
   return gosam_choice
#---#] def ew_gosam_choice:
#---#] def init:

def init():
   """

   We choose which electroweak scheme to follow here. The decision is based
   on what is specified in the model.options line in the input cared
   (for models sm only). These parameters are:

   PARAMETERS:

   1. gosam_choice : the ew_scheme that gosam would choose based on the input
   parameters given in model.option (e.g. mZ=X,mW=Y,alpha=Z => gosam_choice=2)
   2. user_choice  : the ew_scheme the user chose i.e. ewchoose=n
   This is 0 if ewchoose is specified without a number
   3. ewchoose     : A boolean value, True means a choice
   of ew scheme as specified in the file config.f90

   """
   from golem.model.particle import simplify_model
   from golem.model import MODEL_OPTIONS, MODEL_ONES
   global particles, parameters, types, functions

   ew_input = {}

   masses = None
   widths = None

   DEFAULT={}
   DEFAULT['mZ'] = 91.1876
   DEFAULT['mW'] = 80.376
   DEFAULT['alpha'] = 1.0/137.035999679
   DEFAULT['GF'] = 1.16637E-05
   DEFAULT['sw'] = sqrt(0.23120)
   DEFAULT['e'] =  0.3028221202

   keys = MODEL_OPTIONS

   icount = 0
   for key, value in MODEL_OPTIONS.items():
      if key in ["mZ", "mW", "alpha", "GF", "e", "sw"]:
         ew_input[key] = value
         icount += 1

   input_params = set(ew_input.keys())

   ones=golem.model.MODEL_ONES
   eone = False

   gosam_choice = int(ew_gosam_choice(input_params, icount))
   try:
     user_choice = int(keys["users_choice"])
     ewchoose = keys["ewchoose"]
   except:
     user_choice=0
     ewchoose=False
   param=""

   golem.util.tools.debug("GS: %r" % gosam_choice)
   golem.util.tools.debug("user: %r" % user_choice)
   golem.util.tools.debug("ewchoose: %r" % ewchoose)

   for item in DEFAULT.keys():
      if item in ones:
         if item != "e":
            raise Exception("%s is set to one: GoSam cannot handle this EW Scheme" % item)
         else:
            eone = True

   if ewchoose:
      # Substitute parameters from user and fill gaps with defaults once
      for key in ["mZ", "mW", "alpha", "GF", "e", "sw"]:
         if key in MODEL_OPTIONS.keys():
            ew_input[key] = MODEL_OPTIONS[key]
            param+= "%s = %s\n"  % (key, MODEL_OPTIONS[key])
         else:
            ew_input[key] = DEFAULT[key]
      if user_choice != 0:
         warn = "EW scheme was set to ewchoose = %s\n" % user_choice
         warn+= "You specified the following EW parameters:\n"
         if len(param)== 0:
            param = "None!\n"
         warn+= param
         warn+= "We trust you know what you are doing ;-).\n"
         golem.util.tools.warning(warn)
   else:
      if gosam_choice == 0:
         for key in ["mZ","mW","alpha"]:
            ew_input[key] = DEFAULT[key]
      elif gosam_choice == -1:
         # Substitute parameters from user and fill gaps with defaults once
         for key in ["mZ", "mW", "alpha", "GF", "e", "sw"]:
            if key in MODEL_OPTIONS.keys():
               ew_input[key] = MODEL_OPTIONS[key]
               param+= "%s = %s\n"  % (key, MODEL_OPTIONS[key])
            else:
               ew_input[key] = DEFAULT[key]
         keys["ewchoose"] = True
         keys["users_choice"] = '2'
         warn = "EW scheme under-specified.\n"
         warn+= "The number of EW parameters does not allow to select\n"
         warn+= "an EW-scheme. You specified:\n"
         warn+= param
         warn+= "EW parameters are computed from:\n"
         warn+= "mW, mZ, alpha.\n"
         warn+= "This can be changed editing ewchoice in common/config.f90.\n"
         golem.util.tools.warning(warn)
      elif gosam_choice == -2:
         # Substitute parameters from user and fill gaps with defaults once
         for key in ["mZ", "mW", "alpha", "GF", "e", "sw"]:
            if key in MODEL_OPTIONS.keys():
               ew_input[key] = MODEL_OPTIONS[key]
               param+= "%s = %s\n"  % (key, MODEL_OPTIONS[key])
            else:
               ew_input[key] = DEFAULT[key]
         keys["ewchoose"] = True
         keys["users_choice"] = '2'
         warn = "EW scheme over-specified.\n"
         warn+= "The number of EW parameters does not allow to select\n"
         warn+= "an EW-scheme. You specified:\n"
         warn+= param
         warn+= "EW parameters are computed from:\n"
         warn+= "mW, mZ, alpha.\n"
         warn+= "This can be changed editing ewchoice in common/config.f90.\n"
         golem.util.tools.warning(warn)

   for key, value in MODEL_OPTIONS.items():
      if key in parameters:
         try:
            sval = str(value)
            fval = float(sval)
            parameters[key] = sval
         except ValueError:
            golem.util.tools.warning(
            "Model option %s=%r not in allowed range." % (key, value),
            "Option ignored")

      elif key.lower() == "masses":
         if value.strip().lower() == "none":
            masses = []
         else:
            masses = value.split()
      elif key.lower() == "widths":
         if value.strip().lower() == "none":
            widths = []
         else:
            widths = value.split()
   simplify_model(particles, parameters, types, functions, masses, widths)
   init_ew(e_one=eone,**ew_input)

#---#] def init:




#---#[ def fill_dict:

def fill_dict(ad,bd):
   """
   Put this somewhere else
   """
   cd = dict( (key, bd[key]) for key in bd if (key not in ad ))
   for key in cd:
      ad[key] = bd[key]
   return ad
#---#] def fill_dict:
init()
