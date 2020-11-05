#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include "mystic.h"
#include "vroom.h"
#include "uvspecrandom.h"
#include <time.h>
#include "errors.h"
#if HAVE_LIDAR
#include "lidar.h"
#endif
#ifndef PI
#define PI 3.14159265358979323846264338327
#endif
static inline void q42(scadis_struct*q30,int q43);static int q44(pft**
phase_max,int n_phase_max,scadis_struct q30,int q29,int q45,double*mu);static 
inline int q46(double q47,double*mu);static double q48(double q36,double q25,
int q29,int behind_detector,double q49,pft**phase_max,int n_phase_max,double 
q50,double q40,double q51,scadis_struct q30,int q45);int set_vroom_settings(
int vroom,sample_struct*sample,int q7){switch(vroom){case 1:sample->vroom=1;
sample->escape_eps_ddis_upf=0.1;sample->ntupelLE=22;sample->startCP=4;sample->
LEperCP=11;sample->RIS_MS=0;sample->splitter=1;sample->use_p_norm=1;sample->
split_max=3.0;sample->split_min=0.3;sample->n_split_max=5000.0;sample->
n_split_min=0.2;sample->LE_taucrit=3.0;sample->MPdynDDIS=0.1;sample->VIS_CS=1;
sample->vroomreflectalways=0;if(!q7){fprintf(stderr,"\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\124\165\162\156\151\156\147\40\157\156\40\126\122\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\115\56\56\56\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\50\126\141\162\151\141\156\143\145\40\122\145\144\165\143\164\151\157\156\40\117\160\164\151\155\141\154\40\117\160\164\151\157\156\163\40\115\145\164\150\157\144\51\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\111\146\40\171\157\165\40\141\162\145\40\165\163\151\156\147\40\166\162\157\157\155\54\40\160\154\145\141\163\145\40\143\151\164\145\72\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\40\40\122\56\40\102\165\162\141\163\40\141\156\144\40\102\56\40\115\141\171\145\162\40\50\62\60\61\61\51\40\40\40\40\40\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\n"
);fprintf(stderr,"\n");fprintf(stderr,"\52\52\52\40\126\122\117\117\115\40\163\145\164\164\151\156\147\163\72\40\45\144\40\45\144\40\45\144\40\45\144\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\144\40\45\144\n"
,sample->ntupelLE,sample->startCP,sample->LEperCP,sample->splitter,sample->
split_max,sample->split_min,sample->n_split_max,sample->n_split_min,sample->
LE_taucrit,sample->MPdynDDIS,sample->use_p_norm,sample->VIS_CS);}break;case 0:
sample->vroom=0;sample->escape_eps_ddis_upf=0.;sample->ntupelLE=0;sample->
startCP=0;sample->LEperCP=0;sample->RIS_MS=0;sample->splitter=0;sample->
use_p_norm=0;sample->split_max=1e14;sample->split_min=0.;sample->n_split_max=
0.;sample->n_split_min=0.;sample->LE_taucrit=9999.;sample->MPdynDDIS=0.;sample
->VIS_CS=0;sample->vroomreflectalways=0;if(!q7)fprintf(stderr,"\166\162\157\157\155\40\157\146\146\41\40\163\145\164\164\151\156\147\163\72\40\45\144\40\45\144\40\45\144\40\45\144\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\144\40\45\144\n"
,sample->ntupelLE,sample->startCP,sample->LEperCP,sample->splitter,sample->
split_max,sample->split_min,sample->n_split_max,sample->n_split_min,sample->
LE_taucrit,sample->MPdynDDIS,sample->use_p_norm,sample->VIS_CS);break;default:
fprintf(stderr,"\123\124\117\120\41\40\131\157\165\40\141\162\145\40\165\163\151\156\147\40\115\131\123\124\111\103\163\40\154\157\143\141\154\40\145\163\164\151\155\141\164\145\40\164\145\143\150\156\151\161\165\145\54\40\141\156\144\40\171\157\165\162\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\141\162\145\40\163\160\151\153\171\54\n\40\142\165\164\40\171\157\165\40\150\141\166\145\40\156\157\164\40\163\160\145\143\151\146\151\145\144\40\167\150\145\164\150\145\162\40\171\157\165\40\167\141\156\164\40\164\157\40\165\163\145\40\164\150\145\40\166\141\162\151\141\156\143\145\40\162\145\144\165\143\164\151\157\156\40\155\145\164\150\157\144\40\126\122\117\117\115\56\n\40\120\154\145\141\163\145\40\163\160\145\143\151\146\171\40\145\151\164\150\145\162\40\47\155\143\137\166\162\157\157\155\40\157\156\47\40\157\162\40\47\155\143\137\166\162\157\157\155\40\157\146\146\47\40\151\156\40\171\157\165\162\40\151\156\160\165\164\40\146\151\154\145\56\n\40\47\157\156\47\40\151\163\40\162\145\143\157\155\155\145\156\144\145\144\40\146\157\162\40\171\157\165\162\40\143\165\162\162\145\156\164\40\141\160\160\154\151\143\141\164\151\157\156\56\n\105\170\151\164\151\156\147\56\56\56"
);exit(0);}return 0;}int mc_vroom_check_and_verbose(sample_struct*sample,int 
q7,int q8){
#if HAVE_LIDAR
int q52=0;
#endif
if(!(q8)){sample->vroom=0;
#if HAVE_LIDAR
if(sample->LidarLocEst&&sample->LidarLocEst!=q53){q52=set_lidar_settings(
MCLIDAR_NODDIS,sample);if(q52!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\145\164\137\154\151\144\141\162\137\163\145\164\164\151\156\147\163\50\51\n"
,q52);}
#endif
}if(sample->ntupelLE&&!q7)fprintf(stderr,
"\45\144\55\164\165\160\145\154\40\114\105\n",sample->ntupelLE);if(sample->
splitter){if(!(sample->escape_eps_ddis_upf!=0.0||sample->LLE_D_DIS)){fprintf(
stderr,"\127\141\162\156\151\156\147\41\40\143\141\156\47\164\40\165\163\145\40\163\160\154\151\164\164\145\162\40\151\146\40\164\150\145\162\145\40\151\163\40\156\157\40\104\111\123\41\n"
);fprintf(stderr,"\56\56\56\56\56\56\56\56\164\165\162\156\151\156\147\40\157\146\146\40\163\160\154\151\164\164\145\162\56\56\56\56\56\56\56\56\56\56\n"
);sample->splitter=0;}else if(!q7)fprintf(stderr,
"\163\160\154\151\164\164\151\156\147\40\141\142\157\166\145\40\45\145\40\n",
sample->split_max);if(!q7)fprintf(stderr,"\155\141\170\40\163\160\154\151\164\164\151\156\147\40\142\145\154\157\167\40\45\145\40\n"
,sample->n_split_max);}if(sample->LidarLocEst){if(sample->escape){fprintf(
stderr,"\41\41\41\40\123\164\157\160\41\41\41\41\40\131\157\165\40\167\141\156\164\40\155\145\40\164\157\40\145\163\164\151\155\141\164\145\40\145\163\143\141\160\145\40\162\141\144\151\141\156\143\145\163\40\41\41\41\n"
);fprintf(stderr,"\41\41\41\40\141\156\144\40\154\157\143\141\154\40\145\163\164\151\155\141\164\157\162\40\141\164\40\164\150\145\40\163\141\155\145\40\164\151\155\145\41\41\41\40\124\150\151\163\40\151\163\40\40\40\41\41\41\n"
);fprintf(stderr,"\41\41\41\40\156\157\164\40\147\157\151\156\147\40\164\157\40\167\157\162\153\41\41\41\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\41\41\41\n"
);return-1;}}if(!sample->escape&&sample->vroom){fprintf(stderr,"\41\41\41\40\105\162\162\157\162\41\41\41\40\126\122\117\117\115\40\151\163\40\157\156\40\142\165\164\40\105\123\103\101\120\105\40\151\163\40\157\146\146\41\40\n"
);fprintf(stderr,"\41\41\41\40\123\157\155\145\164\150\151\156\147\40\151\163\40\167\162\157\156\147\41\40\103\157\156\164\141\143\164\40\164\150\145\40\144\145\166\145\154\157\160\145\162\163\40\50\143\157\144\145\40\122\102\51\41\n"
);return-1;}if(!q7&&sample->LidarLocEst){fprintf(stderr,"\40\56\56\56\40\162\165\156\156\151\156\147\40\122\125\114\105\123\40\50\114\151\144\141\162\40\105\155\165\154\141\164\157\162\40\45\144\51\40\n"
,sample->LidarLocEst);if(sample->LLE_D_DIS){fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\104\145\164\145\143\164\157\162\40\104\151\162\145\143\164\151\157\156\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\n"
);fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40");if(
sample->LLE_eps_ddis_upf)fprintf(stderr,"\40\165\163\151\156\147\40\120\150\141\163\145\40\106\165\156\143\164\151\157\156\40\50\45\145\51"
,sample->LLE_eps_ddis_upf);if(sample->LLE_eps_ddis_uda)fprintf(stderr,"\40\165\163\151\156\147\40\104\145\164\145\143\164\157\162\40\101\162\145\141\40\50\45\145\51"
,sample->LLE_eps_ddis_uda);if(sample->LLE_eps_fod_dis_phi)fprintf(stderr,"\40\165\163\151\156\147\40\106\151\145\154\144\55\157\146\55\166\151\145\167\40\157\146\40\104\145\164\145\143\164\157\162\40\50\45\145\51"
,sample->LLE_eps_fod_dis_phi);fprintf(stderr,"\n");if(sample->LLE_VIS_FOD)
fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\126\151\162\164\165\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\106\151\145\154\144\55\157\146\55\166\151\145\167\40\117\146\40\104\145\164\145\143\164\157\162\n"
);}if(sample->LLE_VIS_QIDD)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\126\151\162\164\165\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\121\165\141\144\162\141\164\151\143\40\111\156\166\145\162\163\145\40\104\145\164\145\143\164\157\162\40\104\151\163\164\141\156\143\145\n"
);if(sample->LLE_RIS_MAS)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\122\145\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\115\157\154\145\143\165\154\141\162\40\141\156\144\40\101\145\162\157\163\157\154\40\123\143\141\164\164\145\162\151\156\147\n"
);if(sample->LLE_sponti)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\123\160\157\156\164\151\163\160\154\151\164\40\166\145\162\163\151\157\156\40\45\144\n"
,sample->LLE_sponti);
#if HAVE_LIDAR
if(sample->LLE_channels==LIDAR_CHANNEL_RAMAN)fprintf(stderr,"\40\40\40\40\40\101\154\163\157\40\143\141\154\143\165\154\141\164\151\156\147\40\122\141\155\141\156\40\114\151\144\141\162\40\143\150\141\156\156\145\154\56\56\56\n"
);if(sample->LLE_channels==LIDAR_CHANNEL_HSRL)fprintf(stderr,"\40\40\40\40\40\101\154\163\157\40\143\141\154\143\165\154\141\164\151\156\147\40\110\123\122\40\114\151\144\141\162\40\143\150\141\156\156\145\154\56\56\56\n"
);
#endif
if(sample->LLE_turnmax)fprintf(stderr,
"\40\40\40\40\40\124\165\162\156\155\141\170\151\156\147\56\56\56\n");if(
sample->LLE_taumax)fprintf(stderr,
"\40\40\40\40\40\124\141\165\155\141\170\151\156\147\40\45\145\56\56\56\n",
sample->LLE_taumax);if(sample->LE_taucrit)fprintf(stderr,
"	\40\40\40\111\167\141\142\165\143\150\151\156\147\40\45\145\56\56\56\n",
sample->LE_taucrit);fprintf(stderr,"\n");}return 0;}int mc_vroom_prepare(
sample_struct*sample,atmosphere_struct*q9,float q10,int q11,int q12,int q7){
int q54=0,q55=0,q56=0,q57=0,kc=0,q58=0,q59=0,q60=0,q61=0;double q62=0.0,q63=
0.0,q64=0.0;double q65=0.0;int q52=0;int nphamat=1;pft**q66=NULL;int*q67=NULL;
int**q68=NULL;float***q69=NULL,***q70=NULL;double***q71=NULL;double*q72=NULL;
float**q73=NULL,**q74=NULL;double**q75=NULL;int*q76=NULL;int q31=0,q77=0,q78=0
;double q79=0.0;double*F=NULL;double q80=0.0;if(!q7)fprintf(stderr,"\52\52\104\104\111\123\72\40\104\145\146\151\156\151\156\147\40\157\160\164\151\155\141\154\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\146\157\162\40\104\104\111\123\40\56\56\56\n"
);sample->n_phase_max=0;for(q31=1;q31<=q9->n_caoth;q31++){switch(q9->
scatter_type[q31]){case MCSCAT_MOL:case MCSCAT_SKIP:break;case MCSCAT_AER:
sample->n_phase_max+=q11;break;case MCSCAT_HG1:case MCSCAT_HG2:(sample->
n_phase_max)++;break;case MCSCAT_PFT:sample->n_phase_max+=q9->phase[q31]->n;
break;default:fprintf(stderr,"\105\162\162\157\162\54\40\156\157\40\163\165\143\150\40\164\171\160\145\40\45\144\40\157\146\40\163\143\141\164\164\145\162\151\156\147\41\41\41\n"
,q9->scatter_type[q31]);return-1;}}if(sample->n_phase_max!=0){sample->
n_phase_max+=3;q66=calloc(sample->n_phase_max,sizeof(pft*));q67=calloc(sample
->n_phase_max,sizeof(int));q72=calloc(sample->n_phase_max,sizeof(double));q58=
0;for(q31=1;q31<=q9->n_caoth;q31++){switch(q9->scatter_type[q31]){case 
MCSCAT_MOL:q79=0.0;q66[q58]=calloc(1,sizeof(pft));q67[q58]=1;q52=
create_iphase_from_HG(q66[q58++],q79,q7);if(q52!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\162\145\141\164\145\137\151\160\150\141\163\145\137\146\162\157\155\137\110\107\50\51\n"
,q52);if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\155\157\154\145\143\165\154\141\162\40\144\165\155\155\171\72\40\165\163\151\156\147\40\147\75\60\40\141\163\40\151\163\157\164\162\157\160\151\143\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\n"
);break;case MCSCAT_AER:for(q60=0;q60<q11;q60++)if(q9->phase_aer[q60].nphamat
!=0)q66[q58++]=&(q9->phase_aer[q60]);if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\141\154\154\40\117\120\101\103\40\141\145\162\157\163\157\154\40\154\141\171\145\162\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q9->caoth_name[q31]);break;case MCSCAT_HG1:case MCSCAT_HG2:q79=0.0;for(q60=0;
q60<q11;q60++){if(q9->threed[q31][q60]>=1)for(q77=0;q77<q9->Nx;q77++)for(q78=0
;q78<q9->Ny;q78++){if(q79<q9->g1_3D->prof[q31][q60][q77][q78])q79=q9->g1_3D->
prof[q31][q60][q77][q78];if(q79<q9->g2_3D->prof[q31][q60][q77][q78])q79=q9->
g2_3D->prof[q31][q60][q77][q78];}else{if(q79<q9->g1->prof[q31][q60])q79=q9->g1
->prof[q31][q60];if(q79<q9->g2->prof[q31][q60])q79=q9->g2->prof[q31][q60];}}
q66[q58]=calloc(1,sizeof(pft));q67[q58]=1;q52=create_iphase_from_HG(q66[q58++]
,q79,q7);if(q52!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\162\145\141\164\145\137\151\160\150\141\163\145\137\146\162\157\155\137\110\107\50\51\n"
,q52);if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\155\141\170\40\162\137\145\146\146\40\110\107\40\160\150\141\163\145\40\146\143\164\56\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q9->caoth_name[q31]);break;case MCSCAT_PFT:for(q60=0;q60<q9->phase[q31]->n;
q60++)q66[q58++]=q9->phase[q31]->iphase[q60];if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\141\154\154\40\162\137\145\146\146\40\160\150\141\163\145\40\146\143\164\56\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q9->caoth_name[q31]);break;case MCSCAT_SKIP:if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\163\153\151\160\160\151\156\147\40\50\144\165\155\155\171\51\56\56\56\n"
,q9->caoth_name[q31]);break;default:fprintf(stderr,"\105\162\162\157\162\54\40\156\157\40\163\165\143\150\40\164\171\160\145\40\45\144\40\157\146\40\163\143\141\164\164\145\162\151\156\147\41\41\41\n"
,q9->scatter_type[q31]);return-1;}}if(q58>sample->n_phase_max){fprintf(stderr,
"\45\163\45\144\45\163\45\144\45\163","\105\162\162\157\162\41\40\124\150\145\40\163\145\164\40\156\165\155\142\145\162\40\157\146\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\50"
,sample->n_phase_max,"\51\40\163\155\141\154\154\145\162\n\40\164\150\141\156\40\164\150\145\40\156\165\155\142\145\162\40\157\146\40\144\145\146\151\156\145\144\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\50"
,q58,"\51\41\n\40\103\157\156\164\141\143\164\40\122\157\142\145\162\164\40\140\115\145\163\163\171\140\40\102\165\162\141\163\40\146\157\162\40\143\157\155\160\154\141\151\156\164\56\40\105\170\151\164\151\156\147\56\56\56\n"
);return-1;}sample->n_phase_max=q58;if(!q7)fprintf(stderr,"\52\52\104\104\111\123\72\40\146\151\156\141\154\154\171\40\165\163\145\144\40\156\165\155\142\145\162\40\157\146\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\146\157\162\40\104\104\111\123\72\40\45\144\n"
,sample->n_phase_max);q76=calloc(nphamat,sizeof(int));q74=calloc(nphamat,
sizeof(float*));q75=calloc(nphamat,sizeof(double*));q73=calloc(nphamat,sizeof(
float*));q68=calloc(nphamat,sizeof(int*));q69=calloc(nphamat,sizeof(float**));
q71=calloc(nphamat,sizeof(double**));q70=calloc(nphamat,sizeof(float**));for(
q61=0;q61<nphamat;q61++){q68[q61]=calloc(sample->n_phase_max,sizeof(int));q69[
q61]=calloc(sample->n_phase_max,sizeof(float*));q71[q61]=calloc(sample->
n_phase_max,sizeof(double*));q70[q61]=calloc(sample->n_phase_max,sizeof(float*
));for(q58=0;q58<sample->n_phase_max;q58++){q68[q61][q58]=q66[q58]->n[q61];q71
[q61][q58]=q66[q58]->mu[q61];q70[q61][q58]=calloc(q68[q61][q58],sizeof(float))
;for(q60=0;q60<q68[q61][q58];q60++)q70[q61][q58][q60]=(float)q66[q58]->p[
MCSC_MODE_NORMAL][q61][q60];q69[q61][q58]=calloc(q68[q61][q58],sizeof(double))
;for(q60=0;q60<q68[q61][q58];q60++)q69[q61][q58][q60]=acos(q71[q61][q58][q60])
;}}q52=sort_and_add_weighted_phase(sample->n_phase_max,q72,q68,q69,q71,q70,&(
q76),&(q74),&(q75),&(q73),nphamat,1,q7);if(q52!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\157\162\164\137\141\156\144\137\141\144\144\137\167\145\151\147\150\164\145\144\137\160\150\141\163\145\50\51\n"
,q52);sample->phase_max=calloc(1,sizeof(pft*));sample->phase_max[0]=calloc(1,
sizeof(pft));q52=calc_cumulative_table(q75,q73,q76,nphamat,-1.0,sample->
phase_max[0],q10,q7);if(q52!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\143\165\155\165\154\141\164\151\166\145\137\164\141\142\154\145\50\51\n"
,q52);if(q9->nscaDS>1){for(q58=0;q58<sample->n_phase_max;q58++)for(q60=0;q60<
q68[0][q58];q60++)q70[0][q58][q60]=(float)q66[q58]->p[MCSC_MODE_DELTA_SCALE*(
q66[q58]->nscales>1)][0][q60];q52=sort_and_add_weighted_phase(sample->
n_phase_max,q72,q68,q69,q71,q70,&(q76),&(q74),&(q75),&(q73),1,1,q7);if(q52!=0)
return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\157\162\164\137\141\156\144\137\141\144\144\137\167\145\151\147\150\164\145\144\137\160\150\141\163\145\50\51\n"
,q52);for(q59=0;q59<q76[0];q59++)if(q75[0][q59]>q10)break;for(q60=q59+1;q60<
q76[0];q60++)if(q73[0][q60]<q73[0][q59])q73[0][q60]=q73[0][q59];F=calloc(q76[0
],sizeof(double));normalize_phase(q75,q73,F,q76,nphamat,!q7);for(q58=0;q58<
sample->phase_max[0]->n[0];q58++){sample->phase_max[0]->p[
MCSC_MODE_DELTA_SCALE][0][q58]=q73[0][q58];sample->phase_max[0]->F[
MCSC_MODE_DELTA_SCALE][q58]=F[q59];}calc_iphase_coeffs(sample->phase_max[0],
MCSC_MODE_DELTA_SCALE);sample->phase_max[0]->dscale=-999.0;free(F);}for(q58=0;
q58<sample->n_phase_max;q58++)if(q67[q58]){free_iphase(q66[q58]);free(q66[q58]
);}free(q66);free(q67);free(q72);for(q59=0;q59<nphamat;q59++){free(q74[q59]);
free(q75[q59]);free(q73[q59]);}free(q76);free(q74);free(q75);free(q73);for(q61
=0;q61<nphamat;q61++){for(q60=0;q60<sample->n_phase_max;q60++){free(q70[q61][
q60]);free(q69[q61][q60]);}free(q68[q61]);free(q69[q61]);free(q71[q61]);free(
q70[q61]);}free(q68);free(q69);free(q71);free(q70);sample->n_phase_max=1;q80=
0.0;for(q58=0;q58<sample->phase_max[0]->n[0];q58++){if(q80<sample->phase_max[0
]->p[MCSC_MODE_NORMAL][0][q58])q80=sample->phase_max[0]->p[MCSC_MODE_NORMAL][0
][q58];if(q80<1./sample->phase_max[0]->p[MCSC_MODE_NORMAL][0][q58])q80=1./
sample->phase_max[0]->p[MCSC_MODE_NORMAL][0][q58];}if(q80<q2){sample->
n_phase_max=0;free(sample->phase_max[0]);free(sample->phase_max);if(!q7)
fprintf(stderr,"\52\52\116\157\40\156\145\145\144\40\146\157\162\40\104\104\111\123\54\40\164\165\162\156\151\156\147\40\157\146\146\40\126\122\117\117\115\40\50\151\156\40\143\141\163\145\40\151\164\40\167\141\163\40\157\156\51\n"
);}}if(sample->n_phase_max!=0){q52=calloc_hybrid3D_field(&(q9->spiky_box),q9);
if(q52!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\154\157\143\137\150\171\142\162\151\144\63\104\137\146\151\145\154\144\50\51\n"
,q52);for(kc=0;kc<q9->Nz;kc++){if(q9->threed[MCCAOTH_TOT][kc]>=1){for(q31=0;
q31<q9->n_caoth;q31++){if(q9->threed[q31][kc]>=1){for(q54=0;q54<q9->Nx;q54++)
for(q55=0;q55<q9->Ny;q55++)if(q9->ksca3D[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]
->prof[q31][kc][q54][q55]>0.0)q9->spiky_box[kc][q54][q55]=1;}else if(q9->ksca[
MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof[q31][kc]>0.0)for(q54=0;q54<q9->Nx;
q54++)for(q55=0;q55<q9->Ny;q55++)q9->spiky_box[kc][q54][q55]=1;}}else{for(q31=
0;q31<q9->n_caoth;q31++){if(q9->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->
prof[q31][kc]>0.0)q9->spiky_box[kc][0][0]=1;}}}}else{sample->vroom=0;
#if HAVE_LIDAR
if(sample->LidarLocEst&&sample->LidarLocEst!=q53){q52=set_lidar_settings(
MCLIDAR_NODDIS,sample);if(q52!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\145\164\137\154\151\144\141\162\137\163\145\164\164\151\156\147\163\50\51\n"
,q52);}
#endif
}
#ifdef NEWVISQIDD
if(sample->LLE_VIS_QIDD){q62=10.;sample->visqidd_betamax=-q62/sample->lidar[
sample->ili].z_det;sample->visqidd_rmax=-sample->lidar[sample->ili].z_det;
sample->visqidd_facs=sqrt(sample->visqidd_betamax*sample->visqidd_rmax*sample
->visqidd_rmax);fprintf(stderr,"\126\111\123\40\142\145\164\141\155\141\170\40\45\145\40\162\155\141\170\40\45\145\40\146\141\143\40\45\145\n"
,sample->visqidd_betamax,sample->visqidd_rmax,sample->visqidd_facs);
#ifdef NEWRISQIDD
q62=1.;sample->risqidd_betamax=-q62/sample->lidar[sample->ili].z_det;sample->
risqidd_facs=sample->risqidd_betamax*sample->visqidd_rmax*sample->visqidd_rmax
;fprintf(stderr,"\122\111\123\40\142\145\164\141\155\141\170\40\45\145\40\162\155\141\170\40\45\145\40\146\141\143\40\45\145\n"
,sample->risqidd_betamax,sample->visqidd_rmax,sample->risqidd_facs);
#endif
}
#else
if(sample->LLE_VIS_QIDD){q62=0000.0;q63=-q62/sample->lidar[sample->ili].z_det;
for(kc=0;kc<q9->Nz;kc++){q64=(0.5*(q9->Z[kc]+q9->Z[kc+1])-sample->lidar[0].x[2
])/(-sample->lidar[0].dir.dx[2]*sample->lidar[sample->ili].z_det);if(q64>0.0){
if(q64<1.0)q65=q63;else q65=q63/(q64*q64);if(q9->threed[MCCAOTH_TOT][kc]>=1){
for(q54=0;q54<q9->Nx;q54++)for(q55=0;q55<q9->Ny;q55++)for(q56=0;q56<q9->nscaDS
;q56++)for(q57=0;q57<q9->nscaRIS;q57++)if(q65>q9->kext3D[q56][q57][
MCVIS_MODE_QIDD]->prof[MCCAOTH_TOT][kc][q54][q55])q9->kext3D[q56][q57][
MCVIS_MODE_QIDD]->prof[MCCAOTH_TOT][kc][q54][q55]=q65;}else{for(q56=0;q56<q9->
nscaDS;q56++)for(q57=0;q57<q9->nscaRIS;q57++)if(q65>q9->kext[q56][q57][
MCVIS_MODE_QIDD]->prof[MCCAOTH_TOT][kc])q9->kext[q56][q57][MCVIS_MODE_QIDD]->
prof[MCCAOTH_TOT][kc]=q65;}}}}
#endif
return 0;}int mc_vroom_cloning(photon_struct*p,atmosphere_struct*q9,
sample_struct*sample,result_struct*q13,elevation_struct*q14,albedo_struct*q15,
triangular_surface_struct*q16,surftemp_struct*q17,int*q18,int*q19,int q20,int 
q21,float*q22,float*q23,float*q24,int q7){int q52=0;photon_struct*q81=NULL;if(
sample->ntupelLE>1&&!p->isclone&&p->scattercounter>=sample->startCP&&p->
escapescattercounter+sample->LEperCP<p->scattercounter+sample->ntupelLE){q81=
calloc_photon(sample,q9->Nx,q9->Ny,q9->Nz,*q19,q9->nlambda_abs,q9->Nc,q9->
n_caoth);
#ifdef MUCHOUT
if(p->q82==q83)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\143\154\157\156\151\156\147\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode);
#endif
cp_photon_struct(q81,p,sample,q9->n_caoth);q81->isclone=1;q81->wtree=p->wtree;
q81->photon_status=MCSTATUS_SPLIT;q52=photon_journey(q81,q9,sample,q13,q14,q15
,q16,q17,q18,q19,q20,q21,q22,q23,q24,q7,"");if(q52<0){fprintf(stderr,"\105\162\162\157\162\40\45\144\40\157\143\143\165\162\151\156\147\40\151\156\40\160\150\157\164\157\156\137\152\157\165\162\156\145\171\50\51\40\146\157\162\40\160\150\157\164\157\156\40\45\144\40\50\143\154\157\156\145\51\54\40\145\170\151\164\151\156\147\56\56\56\n"
,q52,p->photoncounter);return-1;}destroy_photon(q81,q9->n_caoth);p->
escapescattercounter=p->scattercounter+sample->ntupelLE-1;}return 0;}int 
mc_vroom_splitting_and_rr(photon_struct*p,atmosphere_struct*q9,sample_struct*
sample,result_struct*q13,elevation_struct*q14,albedo_struct*q15,
triangular_surface_struct*q16,surftemp_struct*q17,int*q18,int*q19,int q20,int 
q21,float*q22,float*q23,float*q24,int q7){double q25=0.0,q84=0.0,q85=0.0;int 
q52=0,q86=0,q87=0;
#if HAVE_LIDAR
locest_struct lest;
#endif
photon_struct*q88=NULL;double n_split_max=0.0,n_split_min=0.0;int q89=0;if(
sample->LLE_sponti){if(p->scattercounter<7&&sample->LLE_sponti==1)p->
special_weight*=1.5;else p->special_weight*=1.0+0.1/(p->scattercounter+1);}if(
sample->splitter){if(sample->escape)v_mult_mu(p->dir.dx,sample->rad[0].dir.dx,
&q25);
#if HAVE_LIDAR
if(sample->LidarLocEst){q52=calc_locest_connection(p->x,sample->lidar[sample->
ili].dir.dx,sample->lidar[sample->ili].x,&lest);if(q52<0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\154\157\143\145\163\164\137\143\157\156\156\145\143\164\151\157\156\50\51\n"
,q52);v_mult_mu(p->dir.dx,lest.dir.dx,&q25);}
#endif
q84=get_phase_max(sample->phase_max,sample->n_phase_max,q25,p->DDIS_SC_mode);
if(sample->use_p_norm){q84/=p->p_norm;if(p->p_norm>1.0){if(q84<1.0)p->p_norm*=
q84;if(p->p_norm<1.0)p->p_norm=1.0;}}q85=p->special_weight*p->weight*p->stokes
[0]*exp(-p->tauris)*q84;
#ifdef MUCHOUT
if(p->q82==q83)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\164\145\163\164\40\163\160\154\151\164\164\151\156\147\72\40\162\145\163\164\40\45\145\40\167\145\151\147\150\164\40\45\145\40\160\155\141\170\40\45\145\40\160\156\157\162\155\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,p
->special_weight*p->stokes[0]*exp(-p->tauris),p->weight,q84*p->p_norm,p->
p_norm);
#endif
n_split_max=sample->n_split_max;n_split_min=sample->n_split_min;
#ifdef MUCHOUT
if(p->q82==q83)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\164\145\163\164\40\163\160\154\151\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q85,n_split_max);
#endif
if(q85>sample->split_max&&!(sample->ntupelLE&&!p->isclone&&p->scattercounter>=
sample->startCP)){
#ifdef MUCHOUT
if(p->q82==q83)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\163\160\154\151\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q85,n_split_max);
#endif
q85=q85>n_split_max?n_split_max:q85;q86=(int)q85;p->weight/=(double)q86;q88=
calloc_photon(sample,q9->Nx,q9->Ny,q9->Nz,*q19,q9->nlambda_abs,q9->Nc,q9->
n_caoth);
#ifdef MUCHOUT
if(p->q82==q83)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\163\160\154\151\164\145\163\164\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\144\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q85,q86);
#endif
if(q86>1000){p->spikewarningcounter++;if(p->spikewarningcounter>1){fprintf(
stderr,"\127\141\162\156\151\156\147\41\40\123\160\151\153\145\40\167\141\162\156\151\156\147\40\154\145\166\145\154\40\45\144\40\141\164\40\143\154\157\156\145\40\163\143\141\164\164\145\162\40\157\162\144\145\162\40\45\144\40\163\143\141\164\164\145\162\40\45\144\40\160\150\157\164\157\156\40\45\144\n"
,p->spikewarningcounter,p->clonescattercounter,p->scattercounter,p->
photoncounter);fprintf(stderr,"\40\161\137\163\160\40\45\145\40\167\40\45\145\40\111\60\40\45\145\40\145\170\160\50\55\164\141\165\51\40\45\145\40\120\40\45\145\n"
,p->special_weight,p->weight*(double)q86,p->stokes[0],exp(-p->tauris),q84);}
q89=1;}for(q87=0;q87<q86-1;q87++){cp_photon_struct(q88,p,sample,q9->n_caoth);
q88->wtree=p->wtree;p->photon_status=MCSTATUS_TRAVEL;
#ifdef MUCHOUT
if(p->q82==q83)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\143\157\160\171\156\165\155\142\145\162\72\40\156\137\163\160\154\151\164\40\45\144\40\43\40\45\144\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q86,q87);
#endif
q52=photon_journey(q88,q9,sample,q13,q14,q15,q16,q17,q18,q19,q20,q21,q22,q23,
q24,q7,"");if(q52<0){fprintf(stderr,"\105\162\162\157\162\40\45\144\40\157\143\143\165\162\151\156\147\40\151\156\40\160\150\157\164\157\156\137\152\157\165\162\156\145\171\50\51\40\146\157\162\40\160\150\157\164\157\156\40\45\144\40\50\163\160\154\151\164\40\45\144\51\54\40\145\170\151\164\151\156\147\56\56\56\n"
,q52,p->photoncounter,q87);return-1;}}if(q89){if(p->spikewarningcounter>1)
fprintf(stderr,"\114\145\141\166\151\156\147\40\163\160\151\153\145\40\167\141\162\156\151\156\147\40\154\145\166\145\154\40\45\144\n"
,p->spikewarningcounter);p->spikewarningcounter--;}destroy_photon(q88,q9->
n_caoth);}if(q85<p->weight)q85=p->weight;if(q85<sample->split_min){q85=q85<
n_split_min?n_split_min:q85;if(uvspec_random()>q85)return MCSTATUS_PURGE;p->
weight/=q85;}}return MCSTATUS_DEFAULT;}double get_phase_max(pft**phase_max,int
 n_phase_max,double q25,int SC_mode){int q87=0,q52=0;double q84=0.0,phase=0.0;
for(q87=0;q87<n_phase_max;q87++){q52=get_phase_matrix_pft(phase_max[q87],q25,
SC_mode,1,&phase);if(q52){fct_err_out(q52,"\147\145\164\137\160\150\141\163\145\137\155\141\164\162\151\170\137\160\146\164"
,ERROR_POSITION);return-1.0;}q84+=phase;}return q84/((double)n_phase_max);}
void cp_locest(locest_struct*q26,locest_struct*q19,sample_struct*sample,int 
n_caoth){int q87=0,q31=0,q61=0;cp_direction(&(q26->dir),&(q19->dir));q26->
cosalpha=q19->cosalpha;q26->pdir=q19->pdir;q26->pdir_iso=q19->pdir_iso;if(q26
->pdir_sct!=NULL)for(q31=0;q31<n_caoth+1;q31++)for(q61=0;q61<sample->nstokes;
q61++)q26->pdir_sct[q31][q61]=q19->pdir_sct[q31][q61];q26->dist=q19->dist;q26
->distinv=q19->distinv;q26->r_det=q19->r_det;q26->z_det=q19->z_det;for(q87=0;
q87<3;q87++)q26->x_cc[q87]=q19->x_cc[q87];q26->t_det=q19->t_det;for(q87=0;q87<
3;q87++)q26->x_hit[q87]=q19->x_hit[q87];q26->weight_hit=q19->weight_hit;q26->
in_cone=q19->in_cone;q26->will_hit_cone=q19->will_hit_cone;q26->
behind_detector=q19->behind_detector;q26->will_hit_det_plane=q19->
will_hit_det_plane;q26->lidar_outside_grid=q19->lidar_outside_grid;q26->
hit_det_plane_step=q19->hit_det_plane_step;q26->vis_fod_step=q19->vis_fod_step
;q26->vis_fod_step2=q19->vis_fod_step2;q26->vis_fod_kext=q19->vis_fod_kext;for
(q87=0;q87<3;q87++)q26->hitpoint[q87]=q19->hitpoint[q87];}static inline void 
q42(scadis_struct*q30,int q43){int q87=0;q30->q3=0.0;q30->q4=0.0;q30->q5=2.0;
for(q87=0;q87<3;q87++)q30->dirold_dx[q87]=0.0;q30->d_phi=0.0;q30->epsfac=1.0;
q30->mu_max=1.0;q30->mu_min=-1.0;for(q87=0;q87<q43;q87++)q30->q6[q87]=2.0;for(
q87=0;q87<q43;q87++)q30->F_min[q87]=0.0;}int mc_vroom_prep_DDIS(sample_struct*
sample,photon_struct*p,atmosphere_struct*q9,int*q27,int*q28,int*q29,
locest_struct*lest,scadis_struct*q30){int q87=0;
#if HAVE_LIDAR
int q52=0;double q90=0.0,q91=0.0,q92=0.0;double q93=0.0,q94=0.0;double q95=0.0
;double q96=0.0;
#endif
q42(q30,sample->n_phase_max);if(sample->escape_eps_ddis_upf!=0){if(sample->
ntupelLE&&!p->isclone&&q30->epsfac){if(p->scattercounter<sample->startCP)q30->
epsfac=sample->MPdynDDIS/sample->escape_eps_ddis_upf;else q30->epsfac=0.0;}if(
uvspec_random()<sample->escape_eps_ddis_upf*q30->epsfac)*q27=MCDDIS_UPF;for(
q87=0;q87<3;q87++)lest->dir.dx[q87]=sample->rad[0].dir.dx[q87];}
#if HAVE_LIDAR
if(sample->LLE_D_DIS){q52=calc_locest_connection(p->x,sample->lidar[sample->
ili].dir.dx,sample->lidar[sample->ili].x_cc,lest);if(q52!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\154\157\143\145\163\164\137\143\157\156\156\145\143\164\151\157\156\50\51\n"
,q52);q30->q3=-lest->dist*lest->cosalpha;q30->q4=lest->dist*sqrt(1.-lest->
cosalpha*lest->cosalpha);q90=1./(q30->q3-sample->lidar[sample->ili].z_det);if(
!p->lest.behind_detector){q93=-(q30->q4+sample->lidar[sample->ili].q97)*q90;
q94=acos(lest->cosalpha);q30->q5=cos(atan(q93)-q94);}if(q9->nthreed==0){q95=(p
->x[2]-q9->Z[p->kc])/(q9->Z[p->kc+1]-q9->Z[p->kc]);q30->epsfac*=q95*q9->q98[p
->kc+1]+(1.-q95)*q9->q98[p->kc];if(q30->epsfac<0.)q30->epsfac=0.0;}if(sample->
lidar[sample->ili].cosalpha[0]<lest->cosalpha){if(q30->q3<sample->lidar[sample
->ili].z_det){if(p->lest.behind_detector){fprintf(stderr,"\127\141\162\156\151\156\147\41\40\114\157\147\151\143\141\154\40\163\141\171\163\40\160\150\157\164\157\156\40\151\163\40\142\145\150\151\156\144\40\144\145\164\145\143\164\157\162\54\40\142\165\164\40\151\164\40\151\163\40\151\156\40\146\162\157\156\164\41\n"
);fprintf(stderr,
"\154\157\143\141\164\151\157\156\40\45\145\40\45\145\40\45\145\40\n",p->x[0],
p->x[1],p->x[2]);return-1;}q96=uvspec_random();if(q96<(sample->
LLE_eps_ddis_upf+sample->LLE_eps_ddis_uda)*q30->epsfac)*q27=MCDDIS_UPF;if(q96<
sample->LLE_eps_ddis_uda*q30->epsfac)*q27=MCDDIS_UDA;}else{if(!p->lest.
behind_detector)fprintf(stderr,"\127\141\162\156\151\156\147\41\40\114\157\147\151\143\141\154\40\163\141\171\163\40\160\150\157\164\157\156\40\151\163\40\151\156\40\146\162\157\156\164\40\157\146\40\144\145\164\145\143\164\157\162\54\40\142\165\164\40\151\164\40\151\163\40\142\145\150\151\156\144\41\n"
);*q29=-1;}}else{if(sample->lidar[sample->ili].cosalpha[0]<-lest->cosalpha){*
q29=-1;}else{*q29=1;q96=uvspec_random();if(q96<(sample->LLE_eps_ddis_upf+
sample->LLE_eps_ddis_uda)*q30->epsfac)*q27=MCDDIS_UPF;if(!p->lest.
behind_detector&&q96<sample->LLE_eps_ddis_uda*q30->epsfac)*q27=MCDDIS_UDA;if(*
q27&&uvspec_random()<sample->LLE_eps_fod_dis_phi)*q28=1;q92=-q30->q4/q30->q3;
if(q90>0.)q91=-(q30->q4+sample->lidar[sample->ili].q97)*q90;else q91=-(q30->q4
-sample->lidar[sample->ili].q97)*q90;q95=atan((q91-q92)/(1.0+q91*q92));q30->
mu_max=cos(q95);if(q95<0.)q30->mu_max=-q30->mu_max;if(q95<0.&&q95>-1e-10){q95=
0.0;q30->mu_max=1.0;}q30->mu_min=(q30->q3*sample->lidar[sample->ili].q99-q30->
q4*sample->lidar[sample->ili].q100)*lest->distinv;if(q30->mu_min<-1.0){if(q30
->mu_min<-1.0-q101){fprintf(stderr,"\105\162\162\157\162\41\40\155\165\137\155\151\156\40\75\40\45\145\40\151\156\40\163\143\141\164\164\145\162\151\156\147\50\51\40\151\163\40\165\156\160\150\171\163\151\143\141\154\41\n"
,q30->mu_min);return-1;}q30->mu_min=-1.0;}
#ifdef NOFODDIS
q30->mu_max=1.0;q30->mu_min=-1.0;
#endif
for(q87=0;q87<sample->n_phase_max;q87++){q30->q6[q87]=q102(sample->phase_max[
q87],q30->mu_max,p->DDIS_SC_mode);q30->F_min[q87]=q102(sample->phase_max[q87],
q30->mu_min,p->DDIS_SC_mode);}}}}if(*q29==1){q30->d_phi=sample->lidar[sample->
ili].q100/q30->q4*lest->dist;if(q30->d_phi>1.0)q30->d_phi=1.0;q30->d_phi=q103(
q30->d_phi);}else q30->d_phi=180.;
#endif
return 0;}int mu_scatter_special(atmosphere_struct*q9,photon_struct*p,pft**
phase_max,int n_phase_max,int q31,double*mu,scadis_struct q30,int q27,int q29)
{int q52=0;switch(q27){case MCDDIS_NONE:q52=mu_scatter(q9,p,q31,mu);if(q52!=0)
return err_out("\105\162\162\157\162\54\40\155\165\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q52);break;case MCDDIS_UPF:q52=q44(phase_max,n_phase_max,q30,q29,p->
DDIS_SC_mode,mu);if(q52!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\120\106\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q52);break;case MCDDIS_UDA:q52=q46(q30.q5,mu);if(q52!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\104\101\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q52);break;default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q27);return-1;}return 0;}int random_reflection_special(sample_struct*sample,
atmosphere_struct*q9,elevation_struct*q14,photon_struct*p,pft**phase_max,int 
n_phase_max,double*mu,double*phi,double*q32,scadis_struct*q30,locest_struct 
lest,int q33,int q27,int q28,int q29){int q52=0,q87=0;double q104=0.0;double 
q105=0.0;for(q87=0;q87<3;q87++)q30->dirold_dx[q87]=p->dir.dx[q87];switch(q27){
case MCDDIS_NONE:if(q33==1)random_Lambertian_normal(&(p->dir),q32);else 
random_Isotropic_normal(&(p->dir),q32);v_mult_mu(p->dir.dx,q30->dirold_dx,mu);
break;case MCDDIS_UPF:q52=q44(phase_max,n_phase_max,*q30,q29,p->DDIS_SC_mode,
mu);if(q52!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\120\106\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q52);break;case MCDDIS_UDA:q52=q46(q30->q5,mu);if(q52!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\104\101\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q52);break;default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q27);return-1;}switch(q27){case MCDDIS_NONE:break;case MCDDIS_UPF:case 
MCDDIS_UDA:if(q28)*phi=q30->d_phi*(2.0*uvspec_random()-1.0);else*phi=
sc_Isotropic_phi();if(p->scattercounter==0)q104=p->phi0;else q104=0.0;for(q87=
0;q87<3;q87++)p->dir.dx[q87]=lest.dir.dx[q87];new_direction(*mu,*phi-90.,&(p->
dir),q104);v_mult_mu(p->dir.dx,q32,&q105);if(q105<0.0){p->weight=0.0;}break;
default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q27);return-1;}return 0;}static int q44(pft**phase_max,int n_phase_max,
scadis_struct q30,int q29,int q45,double*mu){int iphase=0;iphase=(int)(
uvspec_random()*((double)n_phase_max-1e-11));if(q30.mu_min==-1.0&&q30.mu_max==
1.0){*mu=sc_mu(phase_max[iphase],1,q45,0.,0.);}else{*mu=sc_mu(phase_max[iphase
],1,q45,q30.q6[iphase],q30.F_min[iphase]);}return 0;}static int inline q46(
double q47,double*mu){double q106=0.0,q107=0.0,q96=0.0;q96=2.*uvspec_random();
q106=q47*q47;q107=2./(1.-q106);if(q96<q107*(q47-q106)){if(q96==0.)*mu=0.;else*
mu=1./(1.+(1.-q47)*(1.-q47)*q107/q96);}else*mu=q106+q96/q107;return 0;}int 
mc_vroom_set_mus_and_phis(sample_struct*sample,photon_struct*p,int q27,
locest_struct lest,scadis_struct q30,double*mu,double*q25,double*q39,double 
phi,double*q40){
#if HAVE_LIDAR
int q52=0;
#endif
switch(q27){case MCDDIS_NONE:v_mult_mu(lest.dir.dx,p->dir.dx,q25);
#if HAVE_LIDAR
if(sample->LLE_VIS_FOD){*q39=sqrt(1.0-*q25**q25);q52=q108(p->dir.dx,lest.dir.
dx,*q39,q40);if(q52)return err_out("\105\122\122\117\122\41\40\144\145\162\151\166\145\137\143\160\150\151\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q52);}
#endif
break;case MCDDIS_UPF:case MCDDIS_UDA:*q25=*mu;v_mult_mu(q30.dirold_dx,p->dir.
dx,mu);if(sample->LLE_VIS_FOD){*q39=sqrt(1.0-*q25**q25);*q40=cosd(phi);}break;
default:fprintf(stderr,"\105\122\122\117\122\41\40\104\104\111\123\163\151\156\147\40\147\151\166\145\163\40\163\157\155\145\164\150\151\156\147\40\163\164\162\141\156\147\145\41\40\45\144\n"
,q27);return-1;}return 0;}int mc_vroom_scattering_calc_phases_and_jacobians(
sample_struct*sample,photon_struct*p,atmosphere_struct*q9,double q34,int q35,
double*q36,double*q37,double*q38){static double**q109=NULL;int q31=0;int q52=0
;if(q35==1){if(q109!=NULL){for(q31=0;q31<=q9->n_caoth;q31++)free(q109[q31]);
free(q109);q109=NULL;}return 0.0;}if(sample->LLE_jacobian||sample->
abs_jacobian)if(q109==NULL){q109=calloc((size_t)q9->n_caoth+1,sizeof(double*))
;for(q31=0;q31<=q9->n_caoth;q31++)q109[q31]=calloc(1,sizeof(double));}q52=
get_phase_matrix_total(q9,p,q34,1,0,sample->spectral_is,sample->
concentration_is,q9->ris_factor,0,q36,q109,q37,&(p->weight));if(q52)return 
fct_err_out(q52,"\147\145\164\137\160\150\141\163\145\137\155\141\164\162\151\170\137\164\157\164\141\154"
,ERROR_POSITION);if(sample->LLE_jacobian)*q38=*q36-q109[MCCAOTH_MOL][0];if(*
q36<=0.0||*q37<=0.0||*q38<0.0){fprintf(stderr,"\105\162\162\157\162\54\40\143\141\154\143\165\154\141\164\151\157\156\40\157\146\40\120\137\156\157\162\155\54\40\120\137\163\160\145\143\54\40\141\156\144\40\120\137\151\163\157\145\156\145\40\144\151\144\40\156\157\164\40\167\157\162\153\41\41\41\n"
);fprintf(stderr,"\120\137\156\157\162\155\40\45\145\40\120\137\163\160\145\143\40\45\145\40\120\137\151\163\157\145\156\145\40\45\145\40\n"
,*q36,*q37,*q38);return-1;}if(sample->LLE_jacobian){if((*q36)!=0.0)for(q31=1;
q31<q9->n_caoth+1;q31++)p->q_jacobian[0][q31-1][p->kc]+=q109[q31][0]/(*q36);if
((*q38)!=0.0)for(q31=1;q31<q9->n_caoth+1;q31++)p->q_jacobian[1][q31-1][p->kc]
+=q109[q31][0]/(*q38);}if(sample->abs_jacobian){for(q31=1;q31<q9->n_caoth+1;
q31++)if(*q36!=0.0){p->q_jacobian[0][q31-1][p->kc]+=q109[q31][0]/(*q36);if(
sample->vroom)p->q_jacobian_sca[q31-1][p->kc]+=q109[q31][0]/(*q36)/get_ksca(q9
,p,q31);}}return 0;}int mc_vroom_DDIS_weight_and_prep_stuff(sample_struct*
sample,atmosphere_struct*q9,double q25,double q39,double q40,int q27,int q29,
double q36,double q37,double q38,locest_struct lest,scadis_struct q30,
photon_struct*p){double cosalpha=0.0;double q110=q37;
#if HAVE_LIDAR
int q52=0;
#endif
if((p->RIS_mode!=MCRIS_MODE_NORMAL)||sample->LLE_D_DIS||sample->
escape_eps_ddis_upf!=0.0||sample->LLE_channels||q9->ris_factor!=1.){if(sample
->escape_eps_ddis_upf!=0.0)q110=q48(q37,q25,0,0,sample->escape_eps_ddis_upf*
q30.epsfac,sample->phase_max,sample->n_phase_max,0,0,0,q30,p->DDIS_SC_mode);
#if HAVE_LIDAR
if(sample->LLE_D_DIS)q110=q48(q37,q25,q29,p->lest.behind_detector,sample->
LLE_eps_ddis_upf*q30.epsfac,sample->phase_max,sample->n_phase_max,sample->
LLE_eps_fod_dis_phi,q40,sample->LLE_eps_ddis_uda*q30.epsfac,q30,p->
DDIS_SC_mode);
#endif
if(!(q110>0.0)){fprintf(stderr,"\105\162\162\157\162\54\40\143\141\154\143\165\154\141\164\151\157\156\40\157\146\40\120\137\163\160\145\143\137\104\40\144\151\144\40\156\157\164\40\167\157\162\153\41\41\41\n"
);fprintf(stderr,"\120\137\163\160\145\143\137\104\40\45\145\40\120\137\163\160\145\143\40\45\145\40\120\137\156\157\162\155\40\45\145\40\n"
,q110,q37,q36);return-1;}p->weight*=q36/q110;
#ifdef MUCHOUT
if(p->q82==q83)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\156\145\167\40\167\145\151\147\150\164\40\45\145\40\146\162\157\155\40\45\145\40\57\40\45\145\40\167\151\164\150\40\155\165\62\40\45\145\40\141\156\144\40\120\163\160\145\143\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,p
->weight,q36,q110,q25,q37);
#endif
p->q_isoene*=q38/q36;}if(sample->LLE_VIS_FOD||sample->LLE_eps_ddis_uda){
v_mult_mu(p->dir.dx,sample->lidar[sample->ili].dir.dx,&cosalpha);p->lest.
hit_det_plane_step=(q30.q3-sample->lidar[sample->ili].z_det)/cosalpha;p->lest.
will_hit_det_plane=(p->lest.hit_det_plane_step>0.0);}
#if HAVE_LIDAR
if(sample->LLE_VIS_FOD){q52=q111(q25,q39,q40,q30.q4,q30.q3,lest.distinv,sample
->lidar[sample->ili].t_det,&p->lest);if(q52!=0)return err_out("\105\162\162\157\162\54\40\143\141\154\143\137\144\151\163\164\141\156\143\145\137\164\157\137\143\157\156\145\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q52);}
#endif
return 0;}static double q48(double q36,double q25,int q29,int behind_detector,
double q49,pft**phase_max,int n_phase_max,double q50,double q40,double q51,
scadis_struct q30,int q45){double q37=0.0;double q112=0.0;double q113=0.0;
double q114=0.0;int q87=0;if(behind_detector){q49+=q51;q51=0.0;}if(q51>0.0)if(
q25>=0.){if(q25>q30.q5)q113=2./(1.-q30.q5*q30.q5);else{q113=(1.-q30.q5)/(1.-
q25);q113*=2./(1.-q30.q5*q30.q5)*q113;}}if(q29==-1)return q36;if(q30.mu_min==-
1.0&&q30.mu_max==1.0){if(q49>0.0)q112=get_phase_max(phase_max,n_phase_max,q25,
q45);q114=1.0;}else{if(q25<=q30.mu_max&&q25>=q30.mu_min)q112=get_phase_max(
phase_max,n_phase_max,q25,q45)*2.0/(q30.q6[q87]-q30.F_min[q87]);if(q40>cosd(
q30.d_phi))q114=1.0-q50*(1.0-180.0/q30.d_phi);else q114=1.0-q50;}q37=(1.0-q49-
q51)*q36+q114*(q49*q112+q51*q113);if(!(q37>0.0)){fprintf(stderr,"\105\162\162\157\162\40\151\156\40\120\137\163\160\145\143\72\40\45\145\40\72\40\145\160\163\137\144\144\151\163\137\165\160\146\40\45\145\40\145\160\163\137\144\144\151\163\137\165\144\141\40\45\145\40\120\137\156\157\162\155\40\45\145\40\146\141\143\137\160\150\151\40\45\145\40\120\137\144\144\151\163\137\165\160\146\40\45\145\40\120\137\144\144\151\163\137\165\144\141\40\45\145\40\155\165\62\40\45\145\40\155\165\155\141\170\40\45\145\40\155\165\155\151\156\40\45\145\n"
,q37,q49,q51,q36,q114,q112,q113,q25,q30.mu_max,q30.mu_min);fprintf(stderr,"\156\137\160\150\141\163\145\137\155\141\170\40\45\144\40\106\137\155\141\170\40\45\145\40\106\137\155\151\156\40\45\145\40\155\165\137\165\144\141\137\144\142\40\45\145\n"
,n_phase_max,q30.q6[0],q30.F_min[0],q30.q5);}
#ifdef MUCHOUT
if(q83==110)fprintf(stderr,"\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\40\143\141\154\143\137\160\137\163\160\145\143\72\40\45\145\40\45\145\40\45\144\40\n"
,q30.mu_min,q30.mu_max,q25<q30.mu_max);
#endif
return q37;}int calloc_hybrid3D_field(int****q41,atmosphere_struct*q9){int kc=
0,q54=0;*q41=calloc((size_t)q9->Nz,sizeof(int**));if(*q41==NULL){fprintf(
stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}for(kc=0;kc<q9->Nz;kc++){if(q9->threed[MCCAOTH_TOT][kc]>=1){(*q41)
[kc]=calloc((size_t)q9->Nx,sizeof(int*));if((*q41)[kc]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}for(q54=0;q54<q9->Nx;q54++){(*q41)[kc][q54]=calloc((size_t)q9->Ny,
sizeof(int));if((*q41)[kc][q54]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}}}else{(*q41)[kc]=calloc((size_t)1,sizeof(int*));if((*q41)[kc]==
NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}q54=0;(*q41)[kc][q54]=calloc((size_t)1,sizeof(int));if((*q41)[kc][
q54]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}}}return 0;}void free_hybrid3D_field(int****q41,atmosphere_struct*
q9){int kc=0,q54=0;for(kc=0;kc<q9->Nz;kc++){if(q9->threed[MCCAOTH_TOT][kc]>=1)
{for(q54=0;q54<q9->Nx;q54++)free((*q41)[kc][q54]);free((*q41)[kc]);}else{free(
(*q41)[kc][0]);free((*q41)[kc]);}}free(*q41);}
