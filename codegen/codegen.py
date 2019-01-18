import sys
#from sage.all import *
from sympy import *
from string import Template

# dict element template
# ('header file','type','typedef')
avail_spaces = { \
	'R':	{'header':'','type':'float','typedef':''}, \
	'R2':	{'header':'Vec.hpp','type':'vec2','typedef':'Vec<float,2>'}}

def gen_cpp_includes_typedefs(main_spaces,min_sl):
	s_includes = ""
	# Includes
	to_include = set([x['header'] for x in min_sl])
	for inc in to_include:
	    if not (inc==''):
	        s_includes+="#include \"headers/{}\"\n".format(inc)
	# Typedefs
	s_typedefs = ""
	typedefs = []
	for x in min_sl:
	    ap = True
	    for y in typedefs:
		if x['typedef']==y['typedef']:
		    ap = False	
	    if ap:
                typedefs.append(x)
	for x in typedefs:
	    if not x['typedef']=='':
		s_typedefs+="typedef {0} {1};\n".format(x['typedef'],x['type'])
	s_typedefs+="\n"
	#define M,Q,TQ
	Mspace = main_spaces[0]
	Qspace = main_spaces[1]
	TQspace = Qspace
	s_typedefs+="typedef {0} M;\ntypedef {1} Q;\ntypedef {2} TQ;".format(Mspace['type'],Qspace['type'],TQspace['type'])
	return {'includes':s_includes,'typedefs':s_typedefs}

def gen_cpp_params(params):
	if len(params)==0:
	    return
	s="struct Params {\n"
	for x in params:
	    s+="\t{0} {1};\n".format(x['type'],x['name'])
	s+="};"
	return {'params':s}

def gen_cpp_disclagsyst(fn,params,meth_type):
	s_step_impl=""
	s_step_impl+="\tM h = js1.base()-js0.base();\n"
	if meth_type=='explicit':
	    Expr = fn
	    Expr = [Expr[i].subs([(p['name'],Symbol('this->m_params.'+p['name'])) for p in params]) for i in range(2)]
	    Expr = [Expr[i].subs([[q,Symbol("js0.pos()")],[v,Symbol("js0.vel()")]]) for i in range(2)]
	    Cexpr = [ccode(Expr[i]) for i in range(2)]
	    s_step_impl+="\tjs1.pos({0});\n\tjs1.vel({1});".format(*Cexpr)
	return {'step_impl':s_step_impl}

def set_spaces_and_types(main_spaces,params):
	global avail_spaces
	l = main_spaces+params
	len_main = len(main_spaces)
	for i in range(len(l)):
	    x = l[i]['space']
	    if not x in avail_spaces.keys():
		sys.exit("Error : unknown space {}\n".format(x))
	    l[i]['type'] = avail_spaces[x]['type']
	min_spacelist = [avail_spaces[y] for y in set([x['space'] for x in l])] #[avail_spaces[x['space']] for x in l]
	return (l[:len_main],l[len_main:],min_spacelist)

def write_cpp_file(filename,ls):
	with open('template.hpp','r') as fi:
	    s = fi.read()
	ts = Template(s)
	res = ts.substitute(**ls)
	with open(filename,'w') as fo:
	    fo.write(res)

base_space = {'space':'R'}
config_space = {'space':'R2'}
meth_type = 'explicit'
n_steps=1000
params=[{'name':'m','value':12.0,'space':'R'}]
main_spaces = [base_space,config_space]
(main_spaces,params,min_spaces) = set_spaces_and_types(main_spaces,params)
ls = {}
ls.update(gen_cpp_includes_typedefs(main_spaces,min_spaces))
ls.update(gen_cpp_params(params))
h,q,v = symbols('h q v')
V=1/q
step=[2+q+h*v,v-h*diff(V,q).subs(q,q+h*v)/symbols('m')]
print(step)
ls.update(gen_cpp_disclagsyst(step,params,meth_type))
print(ls)
filename = 'test.hpp'
write_cpp_file(filename,ls)
