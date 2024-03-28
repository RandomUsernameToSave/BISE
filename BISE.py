from bolos import parser, grid, solver
import numpy as np

class FunctionalFunction(object):
     def __init__(self, func):
             self.func = func
     def __call__(self, *args, **kwargs):
             return self.func(*args, **kwargs)

     def __add__(self, other):
             def summed(*args, **kwargs):
                     return self(*args, **kwargs) + other(*args, **kwargs)
             return FunctionalFunction(summed)

class Species :
    elementSymbol = ''
    fn = FunctionalFunction(lambda x: 0)
    density = 0.
    def __init__(self,elementSymbol, X0):
         self.elementSymbol = elementSymbol
         self.density = X0
         self.fn = FunctionalFunction(lambda x: 0)

    def add_to_function(self,fn):
         self.fn += fn
    
    def get_density(self):
         try:
            return self.density
         except:
              RuntimeError("No initial density passed during creation of element")

    def get_function(self):
         return self.fn

    def get_element(self):
         return self.elementSymbol

class Reaction :
    reactingSpecies = []
    productSpecies = []
    constantReaction = 0.

    def __init__(self, reactingSpecies , productSpecies , constantReaction : float, boltzmanBool= False, boltzmannDat = 'bolsigdb.dat', EN = 120):
        self.reactingSpecies = reactingSpecies
        self.productSpecies = productSpecies
        self.constantReaction = constantReaction
        
        try:
            if boltzmanBool:

                
                gr = grid.LinearGrid(0, 60., 200)
                boltzmann = solver.BoltzmannSolver(gr)
                with open(boltzmannDat) as fp:
                    processes = parser.parse(fp)
                boltzmann.load_collisions(processes)
                if self.reactingSpecies[0].get_element() == 'e':
                    boltzmann.target[self.reactingSpecies[1].get_element()].density = 1.0
                else:
                    boltzmann.target[self.reactingSpecies[0].get_element()].density = 1.0
                boltzmann.kT = 300 * solver.KB / solver.ELECTRONVOLT
                boltzmann.EN = EN * solver.TOWNSEND
                boltzmann.init()
                # Calculate the mean energy according to the first EEDF
                fMaxwell = boltzmann.maxwell(2.0)
                f = boltzmann.converge(fMaxwell, maxn=100, rtol=1e-5)
                mean_energy = boltzmann.mean_energy(f)

                # Set a new grid extending up to 15 times the mean energy.
                # Now we use a quadritic grid instead of a linear one.
                newgrid = grid.QuadraticGrid(0, 15 * mean_energy, 200)

                # Set the new grid and update the internal
                boltzmann.grid = newgrid
                boltzmann.init()

                # Calculate an EEDF in the new grid by interpolating the old one
                finterp = boltzmann.grid.interpolate(f, gr)
                f1 = boltzmann.converge(finterp, maxn=200, rtol=1e-5)
                k = boltzmann.rate(f1, "Ar -> Ar^+")
                self.constantReaction = k
        except:
            raise RuntimeError("Can't calculate boltzmann kinetic constant...")
        
    
    def get_constantReaction(self):
        return self.constantReaction

    def get_reactingSpecies(self):
         return self.reactingSpecies
    
    def get_productSpecies(self):
         return self.productSpecies
    
    def get_function_component(self, list_index):
         
         def fn(allSpeciesOfSolver):
              """Take a vector, return the components associated with the reaction"""
              density = 1
              for ind in list_index:
                   density = density * allSpeciesOfSolver[ind]
              return -self.constantReaction * density
         
         return FunctionalFunction(fn)
    def get_function_component_product(self,list_index):
         def fn(allSpeciesOfSolver):
              """Take a vector, return the components associated with the reaction"""
              density = 1
              
              for ind in list_index:
                   density = density * allSpeciesOfSolver[ind]
              return self.constantReaction * density
         
         return FunctionalFunction(fn)

class Solver : 
    listSpecies = []
    listReaction = []
    initCondition = []

    def get_index_specie(self,specie):
         for i in range(len(self.listSpecies)):
              if self.listSpecies[i] == specie:
                   return i
         else:
              raise Exception("Oups ! specie not integrated into the solver")

    def add_specie(self, Specie):
         self.listSpecies.append(Specie)
    def add_species(self, species):
         self.listSpecies.extend(species)
    
    def add_reaction(self, Reaction):
         self.listReaction.append(Reaction)
    def add_reactions(self, reactions):
         self.listReaction.extend(reactions)
    def summary(self):
         for specie in self.listSpecies:
              print(specie.get_element() + ' ', end="")
         print(' ')
         for react in self.listReaction:
              reacSpecies = react.get_reactingSpecies()
              prodSpecies = react.get_productSpecies()
              k=0
              for rS in reacSpecies:
                   if k== 0:
                        print(specie.get_element() , end="")
                        k = 1
                   else:
                        print(' + ' + specie.get_element() , end="")
              print('  -->  ' , end='')
              k=0
              for rS in reacSpecies:
                   if k== 0:
                        print(specie.get_element() , end="")
                        k= 1
                   else:
                        print(' + ' + specie.get_element() , end="")
              print(' ')

    def init(self):
         for specie in self.listSpecies:
            self.initCondition.append(specie.get_density()) ## This one works well

         for reaction in self.listReaction:
              ## Loop through all reactions
              list_index = []
              for reactSpecie in reaction.get_reactingSpecies():
                   list_index.append(self.get_index_specie(reactSpecie))
              ## ok jusque la
                   
              for reactSpecie in reaction.get_reactingSpecies():  
                   fr = reaction.get_function_component(list_index)
                   reactSpecie.add_to_function(fr)

              for productSpecie in reaction.get_productSpecies():
                   productSpecie.add_to_function(reaction.get_function_component_product(list_index))

         def tot_func(vectorDensity):
              vectorResult = np.zeros_like(np.array(vectorDensity))
              for specie_index in range(len(self.listSpecies)):
                   specie = self.listSpecies[specie_index]
                   fn = specie.get_function()
                   vectorResult[specie_index] = fn(vectorDensity)
             
              return vectorResult
         self.total_function = tot_func
                        
    def solve(self, time_span = [0, 3e-5]):
         from scipy.integrate import solve_ivp, odeint
         t = np.linspace(time_span[0],time_span[1],10)
         y0 = self.initCondition
         def sub_fn(t,v):
            return self.total_function(v)
         #return odeint(sub_fn,y0,t)
         return solve_ivp(sub_fn , time_span , y0)

