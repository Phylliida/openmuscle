import math

# From https://github.com/daniel-haeufle/macroscopic-muscle-model and http://wiki.ifs-tud.de/_media/seminar_3m/3m_2014/haeufle_2014b.pdf

class ContractileElement:
  def __init__(self):
    self.F_max = 1420.0               # F_max in [N] for Extensor (Kistemaker et al., 2006)
    self.l_CEopt =0.092               # optimal length of CE in [m] for Extensor (Kistemaker et al., 2006)
    self.DeltaW_limb_des = 0.35       # width of normalized bell curve in descending branch (Moerl et al., 2012)
    self.DeltaW_limb_asc = 0.35       # width of normalized bell curve in ascending branch (Moerl et al., 2012)
    self.v_CElimb_des = 1.5           # exponent for descending branch (Moerl et al., 2012)
    self.v_CElimb_asc = 3.0           # exponent for ascending branch (Moerl et al., 2012)
    self.A_rel0 = 0.25                # parameter for contraction dynamics: maximum value of A_rel (Guenther, 1997, S. 82)
    self.B_rel0 = 2.25                # parameter for contraction dynmacis: maximum value of B_rel (Guenther, 1997, S. 82)
    # eccentric force-velocity relation:
    self.S_eccentric  = 2.0             # relation between F(v) slopes at v_CE=0 (van Soest & Bobbert, 1993)
    self.F_eccentric  = 1.5           # factor by which the force can exceed F_isom for large eccentric velocities (van Soest & Bobbert, 1993)

class ParallelElasticElement:
  def __init__(self, ce):
    self.L_PEE0   = 0.9                               # rest length of PEE normalized to optimal lenght of CE (Guenther et al., 2007)
    self.l_PEE0   = self.L_PEE0*ce.l_CEopt       # rest length of PEE (Guenther et al., 2007)
    self.v_PEE    = 2.5                               # exponent of F_PEE (Moerl et al., 2012)
    self.F_PEE    = 2.0                               # force of PEE if l_CE is stretched to deltaWlimb_des (Moerl et al., 2012)
    self.K_PEE    = self.F_PEE*( ce.F_max/ ( ce.l_CEopt*(ce.DeltaW_limb_des+1-self.L_PEE0) )**self.v_PEE )
                                                     # factor of non-linearity in F_PEE (Guenther et al., 2007)

class SerialDampingElement:
  def __init__(self, ce):
    self.D_SDE    = 0.3               # xxx dimensionless factor to scale d_SEmax (Moerl et al., 2012)
    self.R_SDE    = 0.01              # minimum value of d_SE normalised to d_SEmax (Moerl et al., 2012)
    self.d_SEmax = self.D_SDE*(ce.F_max*ce.A_rel0)/(ce.l_CEopt*ce.B_rel0)
                                    # maximum value in d_SE in [Ns/m] (Moerl et al., 2012)

class SerialElasticElement:
  def __init__(self):
    self.l_SEE0        = 0.172       # rest length of SEE in [m] (Kistemaker et al., 2006)
    self.DeltaU_SEEnll = 0.0425      # relativ stretch at non-linear linear transition (Moerl et al., 2012)
    self.DeltaU_SEEl   = 0.017       # relativ additional stretch in the linear part providing a force increase of deltaF_SEE0 (Moerl, 2012)
    self.DeltaF_SEE0   = 568.0         # both force at the transition and force increase in the linear part in [N] (~ 40% of the maximal isometric muscle force)

    self.l_SEEnll      = (1 + self.DeltaU_SEEnll)*self.l_SEE0
    self.v_SEE         = self.DeltaU_SEEnll/self.DeltaU_SEEl
    self.KSEEnl        = self.DeltaF_SEE0 / (self.DeltaU_SEEnll*self.l_SEE0)**self.v_SEE
    self.KSEEl         = self.DeltaF_SEE0 / (self.DeltaU_SEEl*self.l_SEE0)

    

    
class Muscle:
  def __init__(self, q):
    self.CE = ContractileElement()
    self.PEE = ParallelElasticElement(self.CE)
    self.SDE = SerialDampingElement(self.CE)
    self.SEE = SerialElasticElement()
    
    self.l_MTC_init = 0.092+0.172  # [m] initial MTC length
    self.q_CE_init  = q          # [] initial muscle activity 0...1
    
    self.l_MTC = self.l_MTC_init
    
    #self.l_CE_init = scipy.optimize.fsolve(lambda x: self.getMuscleForceInit(x, self.l_MTC_init, self.q_CE_init), 0)[0]
    self.l_CE_init = initialConditions.getl_CE_init(q)
    
    
    self.v_CE = 0.0
    self.l_MTC = self.l_MTC_init
    self.l_CE = self.l_CE_init
    self.v_MTC = 0.0 # velocity of muscle, start out at 0
    
    
  
  def step(self, q, dt):
    q = float(q)
    dt = float(dt)
    v_CE = self.compute_v_ce(self.l_CE, self.l_MTC, q)
      
    f_see = self.f_see(self.l_CE, self.l_MTC)
    f_sde = self.f_sde(self.l_CE, q)
    f_mtc = f_see + f_sde
    
    a_mtc = f_mtc / 1000.0
    
    self.v_CE = v_CE
    self.l_CE += self.v_CE*dt
    
    #return str(self.l_MTC) + " " + str(self.v_MTC) + " " + str(self.v_CE) + " " + str(self.l_CE) + " " + str(f_mtc) + " " + str(a_mtc)
    return a_mtc
    
  def f_sde(self, l_CE, q):
    return self.SDE.d_SEmax*((1-self.SDE.R_SDE)*(self.f_ce(l_CE, q) + self.f_pee(l_CE))/self.CE.F_max+self.SDE.R_SDE)*(self.v_MTC-self.v_CE)
  
    
  def f_isom(self, l_CE):
    # Isometric force (Force length relation)
    # Guenther et al. 2007
    if l_CE >= self.CE.l_CEopt: # descending branch:
        F_isom = math.exp( - ( abs( ((l_CE/self.CE.l_CEopt)-1)/self.CE.DeltaW_limb_des ) )**self.CE.v_CElimb_des )
    else: # ascending branch
        F_isom = math.exp( -( abs( ((l_CE/self.CE.l_CEopt)-1)/self.CE.DeltaW_limb_asc ) )**self.CE.v_CElimb_asc )
    return F_isom
    
  def f_pee(self, l_CE):
    # Force of the parallel elastic element
    if l_CE >= self.PEE.l_PEE0:
        F_PEE = self.PEE.K_PEE*(l_CE-self.PEE.l_PEE0)**(self.PEE.v_PEE)
    else: # shorter than slack length
        F_PEE = 0
    return F_PEE
  
  def f_see(self, l_CE, l_MTC):
    # Force of the serial elastic element
    
    l_SEE = l_MTC - l_CE
    if (l_SEE>self.SEE.l_SEE0) and (l_SEE<self.SEE.l_SEEnll): # non-linear part
        F_SEE = self.SEE.KSEEnl*((l_SEE-self.SEE.l_SEE0)**(self.SEE.v_SEE))
    elif l_SEE>=self.SEE.l_SEEnll: # linear part
        F_SEE = self.SEE.DeltaF_SEE0+self.SEE.KSEEl*(l_SEE-self.SEE.l_SEEnll)
    else: # slack length
        F_SEE = 0
    
    return F_SEE
    
    
  
  def f_ce(self, l_CE, q):
    
    aRel, bRel = self.getabRel(l_CE, q, False)
    
    return self.CE.F_max* ( \
      (q*self.f_isom(l_CE)+aRel) / \
      (1 - self.v_CE/(bRel*self.CE.l_CEopt)) \
      - aRel \
      )
  
  
  
  def l_Arel(self, l_CE):
    if l_CE < self.CE.l_CEopt: return 1.0
    else: return self.f_isom(l_CE)
  
  def l_Brel(self):
    return 1.0
    
  def q_Arel(self, q):
    return 1.0/4.0*(1.0+3.0*q)
    
  def q_Brel(self, q):
    return 1.0/7.0*(3.0+4.0*q)
  
  def getabRel(self, l_CE, q, getC):
    aRel = self.CE.A_rel0*self.l_Arel(l_CE)*self.q_Arel(q)
    bRel = self.CE.B_rel0*self.l_Brel()*self.q_Brel(q)
    
    if getC: #self.v_CE > 0:
      f_isom = self.f_isom(l_CE)
      
      f_e = self.CE.F_eccentric
      s_e = self.CE.S_eccentric
      
      aRelC = -f_e*q*f_isom
      bRelC = bRel*(1-f_e)/ \
        (s_e*(1+(aRel/(q*f_isom))))
    
      aRel = aRelC
      bRel = bRelC
    
    return aRel, bRel
  
  
  def compute_v_ce(self, l_CE, l_MTC, q):
  
    f_pee = self.f_pee(l_CE)
    f_isom = self.f_isom(l_CE)
    f_see = self.f_see(l_CE, l_MTC)
    r_se = self.SDE.R_SDE
    f_max = self.CE.F_max
    l_ceOpt = self.CE.l_CEopt
    d_seMax = self.SDE.d_SEmax
    v_mtc = self.v_MTC
    
    aRel, bRel = self.getabRel(l_CE, q, False)
    
    
    d0 = l_ceOpt*bRel*d_seMax*(r_se+(1-r_se)*(q*f_isom+f_pee/f_max))
    c2 = d_seMax*(r_se-(aRel-f_pee/f_max)*(1-r_se))
    c1 = -(c2*v_mtc+d0+f_see-f_pee+f_max*aRel)
    c0 = d0*v_mtc+l_ceOpt*bRel*(f_see-f_pee-f_max*q*f_isom)
    
    v_CE = (-c1-math.sqrt(c1*c1-4*c2*c0))/(2*c2)
    
    if v_CE <= 0:
      aRel, bRel = self.getabRel(l_CE, q, True)
      d0 = l_ceOpt*bRel*d_seMax*(r_se+(1-r_se)*(q*f_isom+f_pee/f_max))
      c2 = d_seMax*(r_se-(aRel-f_pee/f_max)*(1-r_se))
      c1 = -(c2*v_mtc+d0+f_see-f_pee+f_max*aRel)
      c0 = d0*v_mtc+l_ceOpt*bRel*(f_see-f_pee-f_max*q*f_isom)
      
      return (-c1+math.sqrt(c1*c1-4*c2*c0))/(2*c2)
    else:
      return v_CE
    
  
  def getMuscleForceInit(self, l_CE, l_MTC, q):
    F_SEE = self.f_see(l_CE, l_MTC)
    F_CE = q*self.f_isom(l_CE)*self.CE.F_max
    F_PEE = self.f_pee(l_CE)
    
    F_sum = F_SEE-F_CE-F_PEE
    
    return F_sum
        
        
        
        
        
        
        
        
def getl_CE_init(q):
  q = min(max(q, 0.001), 1.0)
  q *= 100
  q = round(q)
  q /= 100
  return initialLengths[q-1]
  
        
        
initialLengths = [
  0.0893969762284,
  0.0891389588719,
  0.0889046126956,
  0.0886891942301,
  0.0884892999683,
  0.0883024013777,
  0.0881265656012,
  0.0879602796499,
  0.0878023352693,
  0.0876517509536,
  0.0875077175724,
  0.0873695594951,
  0.0872367061816,
  0.0871086710171,
  0.0869850352706,
  0.0868654357513,
  0.0867495551782,
  0.0866371145748,
  0.0865278671931,
  0.086421593612,
  0.0863180977442,
  0.086217203557,
  0.0861187523561,
  0.086022600521,
  0.0859286176034,
  0.0858366847203,
  0.0857466931886,
  0.0856585433582,
  0.0855721436094,
  0.0854874094887,
  0.0854042629585,
  0.0853226317445,
  0.0852424487656,
  0.0851636516328,
  0.0850861822091,
  0.0850099862197,
  0.0849350129076,
  0.0848612147265,
  0.0847885470674,
  0.0847169680144,
  0.0846462475631,
  0.0845756343231,
  0.0845050575307,
  0.0844345193694,
  0.0843640220436,
  0.0842935677794,
  0.0842231588239,
  0.0841527974462,
  0.0840824859375,
  0.0840122266114,
  0.0839420218046,
  0.0838718738777,
  0.0838017852163,
  0.0837317582323,
  0.0836617953654,
  0.0835918990854,
  0.0835220718949,
  0.0834523163328,
  0.0833826349788,
  0.0833130304599,
  0.0832435054583,
  0.0831740627232,
  0.0831047050878,
  0.0830354354951,
  0.0829662570414,
  0.0828971730574,
  0.0828281873078,
  0.0827593048977,
  0.0826905305325,
  0.0826218658243,
  0.0825533122074,
  0.0824848711093,
  0.0824165439497,
  0.0823483321402,
  0.0822802370841,
  0.0822122601756,
  0.0821444027999,
  0.0820766663322,
  0.0820090521378,
  0.0819415615714,
  0.081874195977,
  0.0818069566871,
  0.0817398450227,
  0.0816728622929,
  0.0816060097943,
  0.0815392888108,
  0.0814727006134,
  0.0814062464596,
  0.0813399275932,
  0.0812737452442,
  0.0812077006279,
  0.0811417949452,
  0.0810760293822,
  0.0810104051095,
  0.0809449232824,
  0.0808795850406,
  0.0808143915074,
  0.0807493437904,
  0.0806844429803,
  0.0806196901512
]
