# vim: ts=3:sw=3
import ast
import golem
import golem.util.parser
from golem.util.config import Properties

class OLPTemplate(golem.util.parser.Template):
	"""
	Template used to generate files for the Les Houches One-Loop interface.
	"""

	def init_contract(self, contract_file):
		self._contract_file = contract_file
		self._subprocesses = None
		self._subprocesses_main = None
		self._pstack = []
		self._pstack_main = []
		self._listed_flags=set()
		self._listed_flags_length=0


	def init_channels(self, subprocesses, subprocesses_conf):
		self._subprocesses = subprocesses
		self._subprocesses_conf = subprocesses_conf
		
	def init_main_channels(self, subprocesses_main, subprocesses_main_conf):
		self._subprocesses_main = subprocesses_main
		self._subprocesses_main_conf = subprocesses_main_conf

	def subprocesses(self, *args, **opts):
		if "prefix" in opts:
			prefix = opts["prefix"]
		else:
			prefix = ""
		first_name = self._setup_name("first", prefix + "is_first", opts)
		last_name = self._setup_name("last", prefix + "is_last", opts)
		id_name = self._setup_name("id", prefix + "id", opts)
		index_name = self._setup_name("index", prefix + "index",	opts)
		path_name = self._setup_name("path", prefix + "path", opts)
		name_name = self._setup_name("var", prefix + "$_", opts)
		numlegs_name = self._setup_name("num_legs", prefix + "num_legs", opts)
		numhelis_name = self._setup_name("num_helicities",
				prefix + "num_helicities", opts)
		helsum_name = self._setup_name("helsum", prefix + "helsum", opts)
		#inout = self._setup_name("inout",prefix + "inout", opts)
		#print self._contract_file._processes
		#print self._contract_file._proc_res
		#response0= self._contract_file.getProcessResponse(0)
		#self._contract_file.setProcessResponse(response0)

                #for i, inp, out in self._contract_file._processes:
                #    print i, inp, out 
                #    response = self._contract_file.getProcessResponse(0)
                #    print 'response', response

		last = len(self._subprocesses) - 1
		for index, subprocess in enumerate(self._subprocesses):
			props = self._subprocesses_conf[index]
			#print self._contract_file._proc_res[index]
			#i, inp, out= self._contract_file._processes[self._contract_file._proc_res[index][0]]
			#print i, inp, out
			#print index
			props[first_name] = (index == 0)
			props[last_name] = (index == last)
			props[index_name] = index
			props[id_name] = int(subprocess)
			props[name_name] = str(subprocess)
			props[path_name] = subprocess.process_path
			props[numlegs_name] = subprocess.num_legs
			props[numhelis_name] = subprocess.num_helicities
			props[helsum_name] = props.getBooleanProperty('helsum')
			#props[inp] = 

			self._pstack.append(subprocess)
			props.final_extensions=True
			yield props
			self._pstack.pop()
			
			
	def subprocesses_main(self, *args, **opts):
		if "prefix" in opts:
			prefix = opts["prefix"]
		else:
			prefix = ""
		first_name = self._setup_name("first", prefix + "is_first", opts)
		last_name = self._setup_name("last", prefix + "is_last", opts)
		id_name = self._setup_name("id", prefix + "id", opts)
		id_ew_name = self._setup_name("id_ew", prefix + "id_ew", opts)
		index_name = self._setup_name("index", prefix + "index",	opts)
		path_name = self._setup_name("path", prefix + "path", opts)
		name_name = self._setup_name("var", prefix + "$_", opts)
		numlegs_name = self._setup_name("num_legs", prefix + "num_legs", opts)
		numhelis_name = self._setup_name("num_helicities",
				prefix + "num_helicities", opts)
                has_ew = self._setup_name("has_ew",prefix + "has_ew", opts)


                try:
                    last = len(self._subprocesses_main) - 1
                except:
                    last = 0
		for index, subprocess in enumerate(self._subprocesses_main):
			props = self._subprocesses_main_conf[index]
			props[first_name] = (index == 0)
			props[last_name] = (index == last)

			props[index_name] = index
			props[id_name] = int(subprocess)
			props[name_name] = str(subprocess)
			props[path_name] = subprocess.process_path
			props[numlegs_name] = subprocess.num_legs
			props[numhelis_name] = subprocess.num_helicities
			props[has_ew] = subprocess.has_ew

			self._pstack_main.append(subprocess)
			props.final_extensions=True
			yield props
			self._pstack_main.pop()

	def generated_helicities(self, *args, **opts):
		if len(self._pstack) == 0:
			raise TemplateError(
					"[% @for crossings %] outside [% @for subprocesses %]")

		subprocess = self._pstack[-1]
		helis = subprocess.generated_helicities

		last = len(helis) - 1

		if "prefix" in opts:
			prefix = opts["prefix"]
		else:
			prefix = ""

		first_name = self._setup_name("first", prefix + "is_first", opts)
		last_name = self._setup_name("last", prefix + "is_last", opts)
		index_name = self._setup_name("index", prefix + "index",	opts)
		var_name = self._setup_name("var", prefix + "$_", opts)

		if "shift" in opts:
			shift = int(opts[shift])
		else:
			shift = 0

		props = Properties()
		for index, gh in enumerate(helis):
			props[first_name] = (index == 0)
			props[last_name] = (index == last)
			props[index_name] = index 
			props[var_name] = gh + shift

			yield props


	def crossings(self, *args, **opts):
		if "prefix" in opts:
			prefix = opts["prefix"]
		else:
			prefix = ""
		first_name = self._setup_name("first", prefix + "is_first", opts)
		last_name = self._setup_name("last", prefix + "is_last", opts)
		id_name = self._setup_name("id", prefix + "id", opts)
		index_name = self._setup_name("index", prefix + "index",	opts)
		name_name = self._setup_name("var", prefix + "$_", opts)
		channels_name = self._setup_name("channels", prefix + "channels", opts)
		amplitudetype = self._setup_name("amplitudetype",
				prefix + "amplitudetype", opts)
		notreelevel = self._setup_name("notreelevel",
				prefix + "notreelevel", opts)


		include_self = "include-self" in args

		if len(self._pstack) == 0:
			raise TemplateError(
					"[% @for crossings %] outside [% @for subprocesses %]")

		subprocess = self._pstack[-1]


		ids = subprocess.getIDs()
		if not include_self:
			ids.remove(int(subprocess))

		last = len(ids) - 1
		for index, id in enumerate(ids):
			props = Properties()
			props[first_name] = (index == 0)
			props[last_name] = (index == last)
			props[index_name] = index 

			props[id_name] = id
			props[name_name] = subprocess.ids[id]
			if id in subprocess.channels:
				props[channels_name] = subprocess.channels[id]
			else:
				# should happen only if there occured an error before
				props[channels_name]=[]
			props[amplitudetype] = subprocess.getIDConf(id)["olp.amplitudetype"]
			props[notreelevel] = subprocess.getIDConf(id)["olp.notreelevel"]


			yield props
			
	def crossings_main(self, *args, **opts):
		if "prefix" in opts:
			prefix = opts["prefix"]
		else:
			prefix = ""
		first_name = self._setup_name("first", prefix + "is_first", opts)
		last_name = self._setup_name("last", prefix + "is_last", opts)
		id_name = self._setup_name("id", prefix + "id", opts)
		id_ew_name = self._setup_name("id_ew", prefix + "id_ew", opts)
		index_name = self._setup_name("index", prefix + "index",	opts)
		name_name = self._setup_name("var", prefix + "$_", opts)
		channels_name = self._setup_name("channels", prefix + "channels", opts)
		amplitudetype = self._setup_name("amplitudetype",
				prefix + "amplitudetype", opts)
		notreelevel = self._setup_name("notreelevel",
				prefix + "notreelevel", opts)
                inout = self._setup_name("inout",prefix + "inout", opts)
                response = self._setup_name("response",prefix +"response", opts)


		include_self = "include-self" in args

		if len(self._pstack_main) == 0:
			raise TemplateError(
					"[% @for crossings %] outside [% @for subprocesses %]")

		subprocess = self._pstack_main[-1]

		ids = subprocess.getIDs()
		if not include_self:
			ids.remove(int(subprocess))

                processes={}

                id_list=None
                j=-1
                while self._subprocesses_main_conf[j]["id_list"]==None:
                  j=j-1
                id_list = ast.literal_eval(self._subprocesses_main_conf[j]["id_list"])                  

                for i, inp, out in self._contract_file.processes():
                    processes[id_list[id_list.keys()[i]][0]]=(inp,out)                

		last = len(ids) - 1
		for index, id in enumerate(ids):
			props = Properties()
			props[first_name] = (index == 0)
			props[last_name] = (index == last)
			props[index_name] = index 
			props[id_name] = id
			props[name_name] = subprocess.ids[id]
			if id in subprocess.channels:
				props[channels_name] = subprocess.channels[id]
				res = self._contract_file.getProcessResponse(id)
				inp, out= processes[res[0]]
				props[inout]="%s -> %s" % (" ".join(map(str,inp))," ".join(map(str,out)))
				props[response]=self._contract_file.getProcessResponse(id)
			else:
				# should happen only if there occured an error before
				props[channels_name]=[]
			props[amplitudetype] = subprocess.getIDConf(id)["olp.amplitudetype"]
			props[notreelevel] = subprocess.getIDConf(id)["olp.notreelevel"]


			yield props

	def flags_filter(self, *args, **opts):
		assert args == ("$_",)
		flags=self._stack[-1]["$_"]
		ret=""
		for f in flags.split():
			if f not in self._listed_flags:
				self._listed_flags_length += len(f)
				if self._listed_flags_length > 80:
					self._listed_flags_length=0
					ret = ret + "\\\n "
					self._listed_flags_length += len(f)
				ret = ret + f + " "
				self._listed_flags.add(f)
		if self._stack[-1]["is_last"]=="True":
			ret=ret.rstrip()
		return ret

	def reset_flags_filter(self, *args, **opts):
		self._listed_flags=set()
		self._listed_flags_length=0
		return ""
