using System;

// From https://github.com/daniel-haeufle/macroscopic-muscle-model and http://wiki.ifs-tud.de/_media/seminar_3m/3m_2014/haeufle_2014b.pdf

public class Muscle
{
        ContractileElement CE;
        ParallelElasticElement PEE;
        SerialDampingElement SDE;
        SerialElasticElement SEE;


        double l_MTC;
        double v_MTC;

        double l_CE;
        double v_CE;


        public Muscle(double q)
        {
            CE = new ContractileElement();
            PEE = new ParallelElasticElement(CE);
            SDE = new SerialDampingElement(CE);
            SEE = new SerialElasticElement();


            double l_MTC_init = 0.092 + 0.172;  // [m] initial MTC length
            double q_CE_init = q;          // [] initial muscle activity 0...1
            double l_CE_init = getl_CE_init(q_CE_init);

            l_MTC = l_MTC_init;
            l_CE = l_CE_init;

            v_MTC = 0.0;
            v_CE = 0.0;
        }

        public double step(double q, double dt)
        {
            q = (float)(q);
            dt = (float)(dt);
            double v_CE = this.compute_v_ce(this.l_CE, this.l_MTC, q);

            double f_see = this.f_see(this.l_CE, this.l_MTC);
            double f_sde = this.f_sde(this.l_CE, q);
            double f_mtc = f_see + f_sde;

            double a_mtc = f_mtc / 1000.0;


            this.v_CE = v_CE;
            this.l_CE += this.v_CE * dt;

            return a_mtc;
        }

        double f_sde(double l_CE, double q)
        {
            return this.SDE.d_SEmax * ((1 - this.SDE.R_SDE) * (this.f_ce(l_CE, q) + this.f_pee(l_CE)) / this.CE.F_max + this.SDE.R_SDE) * (this.v_MTC - this.v_CE);
        }

        double f_isom(double l_CE)
        {
            // Isometric force (Force length relation)
            // Guenther et al. 2007
            double F_isom;
            if (l_CE >= this.CE.l_CEopt) // descending branch:
                F_isom = Math.Exp(-Math.Pow((Math.Abs(((l_CE / this.CE.l_CEopt) - 1) / this.CE.DeltaW_limb_des)), this.CE.v_CElimb_des));
            else // ascending branch
                F_isom = Math.Exp(-Math.Pow((Math.Abs(((l_CE / this.CE.l_CEopt) - 1) / this.CE.DeltaW_limb_asc)), this.CE.v_CElimb_asc));
            return F_isom;
         }

        double f_pee(double l_CE)
        {
            // Force of the parallel elastic element
            double F_PEE;

            if (l_CE >= this.PEE.l_PEE0)
                F_PEE = this.PEE.K_PEE * Math.Pow((l_CE - this.PEE.l_PEE0), (this.PEE.v_PEE));
            else // shorter than slack length
                F_PEE = 0;
            return F_PEE;
         }

        double f_see(double l_CE, double l_MTC)
        {
            // Force of the serial elastic element

            double l_SEE = l_MTC - l_CE;
            double F_SEE;
            if (l_SEE > this.SEE.l_SEE0 && l_SEE < this.SEE.l_SEEnll) // non-linear part
                F_SEE = this.SEE.KSEEnl * (Math.Pow((l_SEE - this.SEE.l_SEE0), (this.SEE.v_SEE)));
            else if (l_SEE>= this.SEE.l_SEEnll) // linear part
                F_SEE = this.SEE.DeltaF_SEE0 + this.SEE.KSEEl * (l_SEE - this.SEE.l_SEEnll);
            else // slack length
                F_SEE = 0;

            return F_SEE;
         }

        double f_ce(double l_CE, double q)
        {


            double aRel, bRel;
            this.getabRel(l_CE, q, false, out aRel, out bRel);


            return this.CE.F_max * (
              (q * this.f_isom(l_CE) + aRel) /
              (1 - this.v_CE / (bRel * this.CE.l_CEopt))
              - aRel
              );
         }

        double l_Arel(double l_CE)
        {
            if (l_CE < this.CE.l_CEopt) return 1.0;
            else return this.f_isom(l_CE);
        }

        double l_Brel()
        {
            return 1.0;
        }

        double q_Arel(double q)
        {
            return 1.0 / 4.0 * (1.0 + 3.0 * q);
        }

        double q_Brel(double q)
        {
            return 1.0 / 7.0 * (3.0 + 4.0 * q);
        }

        void getabRel(double l_CE, double q, bool getC, out double aRel, out double bRel)
        {
            aRel = this.CE.A_rel0 * this.l_Arel(l_CE) * this.q_Arel(q);
            bRel = this.CE.B_rel0 * this.l_Brel() * this.q_Brel(q);


            if (getC) //this.v_CE > 0:
            {
                double f_isom = this.f_isom(l_CE);


                double f_e = this.CE.F_eccentric;
                double s_e = this.CE.S_eccentric;


                double aRelC = -f_e * q * f_isom;
                double bRelC = bRel * (1 - f_e) /
                   (s_e * (1 + (aRel / (q * f_isom))));


                aRel = aRelC;
                bRel = bRelC;
            }
         }

        double compute_v_ce(double l_CE, double l_MTC, double q)
        {


            double f_pee = this.f_pee(l_CE);
            double f_isom = this.f_isom(l_CE);
            double f_see = this.f_see(l_CE, l_MTC);
            double r_se = this.SDE.R_SDE;
            double f_max = this.CE.F_max;
            double l_ceOpt = this.CE.l_CEopt;
            double d_seMax = this.SDE.d_SEmax;
            double v_mtc = this.v_MTC;
            double aRel, bRel;
            this.getabRel(l_CE, q, false, out aRel, out bRel);



            double d0 = l_ceOpt * bRel * d_seMax * (r_se + (1 - r_se) * (q * f_isom + f_pee / f_max));
            double c2 = d_seMax * (r_se - (aRel - f_pee / f_max) * (1 - r_se));
            double c1 = -(c2 * v_mtc + d0 + f_see - f_pee + f_max * aRel);
            double c0 = d0 * v_mtc + l_ceOpt * bRel * (f_see - f_pee - f_max * q * f_isom);


            double v_CE = (-c1 - Math.Sqrt(c1 * c1 - 4 * c2 * c0)) / (2 * c2);


            if (v_CE <= 0)
            {
                this.getabRel(l_CE, q, true, out aRel, out bRel);
                d0 = l_ceOpt * bRel * d_seMax * (r_se + (1 - r_se) * (q * f_isom + f_pee / f_max));
                c2 = d_seMax * (r_se - (aRel - f_pee / f_max) * (1 - r_se));
                c1 = -(c2 * v_mtc + d0 + f_see - f_pee + f_max * aRel);
                c0 = d0 * v_mtc + l_ceOpt * bRel * (f_see - f_pee - f_max * q * f_isom);


                return (-c1 + Math.Sqrt(c1 * c1 - 4 * c2 * c0)) / (2 * c2);
            }
            else
            {
                return v_CE;
            }
        }


        double getMuscleForceInit(double l_CE, double l_MTC, double q)
        {
            double F_SEE = this.f_see(l_CE, l_MTC);
            double F_CE = q * this.f_isom(l_CE) * this.CE.F_max;
            double F_PEE = this.f_pee(l_CE);


            double F_sum = F_SEE - F_CE - F_PEE;


            return F_sum;
      }


    private double getl_CE_init(double q)
    {
        q = Math.Min(Math.Max(q, 0.001), 1.0);
        q *= 100;
        q = Math.Round(q);
        return initialLengths[(int)(q - 1)];
    }



   static double[] initialLengths = new double[] {
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
   };
}


class ContractileElement
{
    public double F_max = 1420.0;               // F_max in [N] for Extensor (Kistemaker et al., 2006)
    public double l_CEopt = 0.092;              // optimal length of CE in [m] for Extensor (Kistemaker et al., 2006)
    public double DeltaW_limb_des = 0.35;       // width of normalized bell curve in descending branch (Moerl et al., 2012)
    public double DeltaW_limb_asc = 0.35;       // width of normalized bell curve in ascending branch (Moerl et al., 2012)
    public double v_CElimb_des = 1.5;           // exponent for descending branch (Moerl et al., 2012)
    public double v_CElimb_asc = 3.0;           // exponent for ascending branch (Moerl et al., 2012)
    public double A_rel0 = 0.25;                // parameter for contraction dynamics: maximum value of A_rel (Guenther, 1997, S. 82)
    public double B_rel0 = 2.25;                // parameter for contraction dynmacis: maximum value of B_rel (Guenther, 1997, S. 82)

    // eccentric force-velocity relation:
    public double S_eccentric = 2.0;            // relation between F(v) slopes at v_CE=0 (van Soest & Bobbert, 1993)
    public double F_eccentric = 1.5;            // factor by which the force can exceed F_isom for large eccentric velocities (van Soest & Bobbert, 1993)
}

class ParallelElasticElement
{

    public ParallelElasticElement(ContractileElement ce)
    {
        l_PEE0 = L_PEE0 * ce.l_CEopt;
        K_PEE = F_PEE * (ce.F_max / Math.Pow((ce.l_CEopt * (ce.DeltaW_limb_des + 1 - L_PEE0)), v_PEE));
    }


    public double L_PEE0 = 0.9;                               // rest length of PEE normalized to optimal lenght of CE (Guenther et al., 2007)
    public double l_PEE0;                                     // rest length of PEE (Guenther et al., 2007)
    public double v_PEE = 2.5;                                // exponent of F_PEE (Moerl et al., 2012)
    public double F_PEE = 2.0;                                // force of PEE if l_CE is stretched to deltaWlimb_des (Moerl et al., 2012)
    public double K_PEE;                                      // factor of non-linearity in F_PEE (Guenther et al., 2007)
}

class SerialDampingElement
{
    public SerialDampingElement(ContractileElement ce)
    {
        d_SEmax = D_SDE * (ce.F_max * ce.A_rel0) / (ce.l_CEopt * ce.B_rel0);
    }

    public double d_SEmax;                   // maximum value in d_SE in [Ns/m] (Moerl et al., 2012)
    public double D_SDE = 0.3;               // xxx dimensionless factor to scale d_SEmax (Moerl et al., 2012)
    public double R_SDE = 0.01;              // minimum value of d_SE normalised to d_SEmax (Moerl et al., 2012)
}


class SerialElasticElement
{
    public SerialElasticElement()
    {
        l_SEEnll = (1 + DeltaU_SEEnll) * l_SEE0;
        v_SEE = DeltaU_SEEnll / DeltaU_SEEl;
        KSEEnl = DeltaF_SEE0 / Math.Pow((DeltaU_SEEnll * l_SEE0), v_SEE);
        KSEEl = DeltaF_SEE0 / (DeltaU_SEEl * l_SEE0);
    }

    public double l_SEE0 = 0.172;              // rest length of SEE in [m] (Kistemaker et al., 2006)
    public double DeltaU_SEEnll = 0.0425;      // relativ stretch at non-linear linear transition (Moerl et al., 2012)
    public double DeltaU_SEEl = 0.017;         // relativ additional stretch in the linear part providing a force increase of deltaF_SEE0 (Moerl, 2012)
    public double DeltaF_SEE0 = 568.0;         // both force at the transition and force increase in the linear part in [N] (~ 40% of the maximal isometric muscle force)

    public double l_SEEnll;
    public double v_SEE;
    public double KSEEnl;
    public double KSEEl;
}

