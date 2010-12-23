/**********************************************************************/
void dfint(Int_t narg, Double_t *arg, Double_t *nent, Double_t *ent, Double_t *table)
{
    Double_t kd = 1;
    Double_t m = 1;
    Int_t ncomb[narg];
    Int_t ja = 0;
    Int_t jb, k;
    Int_t i = 0;
    while (i <= narg) {
        ncomb[i] = 1;
        jb = ja + nent[i]; // VERIFICAR SI LOS LIMITES ESTAN CORRECTOS!!!!
        Int_t j = ja;
        while (j <= jb) {
            if (arg[i] <= ent[j]) {
                if (j != ja) {
                    jr = j - 1;
                    d[i] = (ent[j] - arg[i]) / (ent[j] - ent[jr]);
                    ient[i] = j - ja;
                    kd = kd + ient[i] * m;
                    m = m * nent[i];
                    ja = jb + 1;
                    break;
                } else {
                    j++;
                    jr = j - 1;
                    d[i] = (ent[j] - arg[i]) / (ent[j] - ent[jr]);
                    ient[i] = j - ja;
                    kd = kd + ient[i] * m;
                    m = m * nent[i];
                    ja = jb + 1;
                    break;
                }

            }
            j++;
        }
        j = jb;
        j++;
        d[i] = (ent[j] - arg[i]) / (ent[j] - ent[jr]);
        ient[i] = j - ja;
        kd = kd + ient[i] * m;
        m = m * nent[i];
        ja = jb + 1;
        i++;
    }
    dfint = 0.;
    while (1) {
        fac = 1.;
        iadr = kd;
        ifadr = 1.;
        //cambiado a 0 dado el indice de los arreglos en fortran
        i = 0;
        while (i <= narg) {
            if (ncomb[i] == 0) {
                fac = fac * d[i];
                iadr  = iadr - ifadr;
                ifadr = ifadr * nent[i];

            } else {
                fac = fac * (1. -  d[i]);
                ifadr = ifadr * nent[i];
            }
            i++;
        }
        dfint = dfint + fac * table[adr];
        il = narg;

        while (ncomb[il] == 0) {
            il--;
            if (il == 0) return;
        }

        ncomb[il] = 0;
        if (il == narg) {
            continue;
        } else {
            il++;
            k = il;
            while (k <= narg) {
                ncomb[k] = 1;
                k++;
            }
            continue;
        }
    }

}
