
/* This code is adapted and modified from
Rosenbaum, R., Smith, M. A., Kohn, A., Rubin, J. E., & Doiron, B. (2017). The spatial structure of correlated neuronal variability. Nature neuroscience, 20(1), 107.
https://github.com/RobertRosenbaum/SpatialSpikingNeuralNetworks/tree/master
*/

/* Inputs:
   sx is the feedforward presynaptic spike trains, which should be 3xNsx where Nsx is the number of spikes.
      sx(1,:) is spike times, sx(2,:) is the x-index of the neuron and sx(3,:) is the y-index
   Nx1 is the number of neurons in the feedforward layer in each direction.
      sx(2,:) and sx(3,:) should be integers between 1 and Nx1 inclusively.
   Ne1 and Ni1 are the numbers of exc and inh neurons in each direction
      So there should be Nx==Nx1^2 neurons in the ffwd layer in all,
      Ne==Ne1^2 exc neurons in the recurrent network and Ni==Ni1^2 inh neurons
   Jab is the synaptic strength of connections from b=e,i,x to a=e,i.
   Kab is the number of projections from each cell in pop b=e,i,x to all cells in pop a=e,i
   betaab is the "width" of connections from a to b (i.e., the std of the gaussian)
      It is in units of neuron indices, so it should be proportional to Na1
   C,gl,vl,DeltaT,VT,tref,Vth,Vre,Vlb are LIF neuron params
      They are each 2x1 vectors, for exc and inh neurons separately.
      For example, C(1) is the capacitance of exc neurons and C(2) of inh neurons
   tausynb is the time-constant of the synapses from neurons in population b=x,e,i.
     post-synaptic currents are of the form (1/tausynb)*exp(-t/tausynb) where t>0 is
     the time evolved since the presynaptic spike.
   V0 is vector of membrane potential initial conditions.
      It should be Nx1 where N=Ne+Ni is the number of neurons in the recurrent network
      The first Ne elements are for exc neurons, the last Ni for inh neurons
   T is the total simulation time
   dt is time bin size
   maxns is maximum number of spikes allowed for all neurons together
   Irecord is a 2x(Nrecord) matrix indicating for which neurons we should record
     the synaptic inputs and membrane potential.
     Irecord(1,:) is the index in the x-direction and Irecord(2,:) in the y-direction.
     Excitatory neurons are the first Ne1 in each direction and inhibitory neurons the
     next Ni1.  For example, it Ne1=100 then Irecord(1,j)=150, Irecord(2,j)=130 means
     that the jth neuron recorded will be the inhibtiory neuron at coordinates
     (50,30)

   Outputs:
   s is a 3x(maxns) matrix of spikes
      s(1,:) contains spike times.
      s(2,:) contains indices of neurons that spike in the x-direction
      s(3,:) contains indices of neurons that spike in the y-direction
      Inhibitory neurons are given negative indices, exc neurons positive.
      For example the jth spike occurs at time s(1,j).
      If s(2,j)==-50 and s(3,j)==-10 then it is an inhibitory neuron spike
      at spatial coordinates (50, 10) or (50/Ni1, 10/Ni1) on the unit square.
      When there are fewer than maxns spikes, extra space in s will be filled
      with zeros.  This should be truncated in matlab by writing s=s(:,s(1,:)>0);
   Ix,Ie,Ii,v are the recorded
   feedforward synaptic input, exc synaptic input, inh synaptic input and voltage
   respectively.
 */

#include "mex.h"
#include "math.h"
#include "stdlib.h"
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef srand48
#define srand48(s) srand(s)
#endif

#ifndef drand48
#define drand48() (((double)rand()) / ((double)RAND_MAX))
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int printfperiod, kk, Ntref[2], Ni, Ni1, Ne2, Ni2, Kee, Kie, Kei, Kii, Kex, Kix, Ne, Ne1, j, j1, j2, k, i, N, Nt, m1, m2, maxns, ns, flagei, Nrecord, Nrecord_I, jj, nskiprecord, flag, tempflag, Nsx, Nx1, Nx, D;
    double dt, *s, *v, *v0, *JnextE, *JnextI, *JnextX, *alphax, *alphae, *alphai, tausyne, tausyni, *alphaxr, *alphaer, *alphair, *vr, *sx, Jex, Jix, betaex, betaix, tausynx, *total_ffwd_input, *total_recurrent_ex_input, *total_recurrent_inh_input, *total_input, *input_count, *input, *JnextE1, *JnextI1, *JnextX1, *alphax1, *alphae1, *alphai1, *alphaxr1, *alphaer1, *alphair1, *JnextE2, *JnextI2, *JnextX2, *alphax2, *alphae2, *alphai2, *alphaxr2, *alphaer2, *alphair2, *JnextE3, *JnextI3, *JnextX3, *alphax3, *alphae3, *alphai3, *alphaxr3, *alphaer3, *alphair3, *JnextE4, *JnextI4, *JnextX4, *alphax4, *alphae4, *alphai4, *alphaxr4, *alphaer4, *alphair4, *input_NL, *input_dendrite_1_before, *input_dendrite_1_after, *input_NL_after;
    double Jee, Jei, Jie, Jii, betaee, betaei, betaie, betaii, T, *Irecord, *Irecord_I, *C, *Vleak, *DeltaT, *VT, *tref, *gl, *Vth, *Vre, *Vlb, xloc, yloc;
    int *tempWee, *tempWei, *tempWie, *tempWii, *tempWex, *tempWix, *Wee1, *Wee2, *Wei1, *Wei2, *Wie1, *Wie2, *Wii1, *Wii2, *Wex1, *Wex2, *Wix1, *Wix2, *refstate, iXspike, jspike, postcell, *Dee, *Dei, *Die, *Dii, *Dex, *Dix, *x_coords_e, *y_coords_e, *x_coords_i, *y_coords_i, *x_coords_x, *y_coords_x;
    mxArray *temp1, *temp2, *temp0, *temp3, *temp4, *temp5, *temp6, *temp7, *temp8, *temp9, *temp10, *temp11, *temp12, *temp13, *temp14, *temp15, *temp16, *temp17, *temp18, *temp19, *temp20, *temp21, *temp22;

    double dendrite_length = 1.0;
    D = 4;
    double lambda = 0.5; // 0.7;
    double centered_sigmoid = 0.0;
    double input_d1 = 0.0;

    // double **dendriteArray = mxMalloc(N * sizeof(double *));
    // for (int i = 0; i < N; i++)
    // {
    //     dendriteArray[i] = mxMalloc(D * sizeof(double));
    // }

    /* Seed random number generator */
    /* Change this if you want */
    /* You could alternatively pass in a seed */
    srand48(10);

    void CircRandNfun(int *, double, double, int, int, int);
    void CircRandNfunTest(int *, int *, double, double, int, int, int);
    int closest_dendrite(double, double, double, double, int, int);
    double applySigmoid(double, double, double);
    void save1DArrayToCSV(const char *, int *, int);

    /******
     * Import variables from matlab
     * This is messy looking and is specific to mex.
     * Ignore if you're implementing this outside of mex.
     *******/

    sx = mxGetPr(prhs[0]);
    m1 = mxGetM(prhs[0]);
    Nsx = mxGetN(prhs[0]);
    if (m1 != 3)
    {
        mexErrMsgTxt("sx should be Nsxx3");
    }

    /* Number of neurons in the ffwd layer in each direction. */
    Nx1 = (int)mxGetScalar(prhs[1]);

    /* Number of exc neurons in each direction. */
    Ne1 = (int)mxGetScalar(prhs[2]);

    /* Number of inh neurons in each direction. */
    Ni1 = (int)mxGetScalar(prhs[3]);

    Jex = mxGetScalar(prhs[4]);
    Jix = mxGetScalar(prhs[5]);

    Jee = mxGetScalar(prhs[6]);
    Jei = mxGetScalar(prhs[7]);
    Jie = mxGetScalar(prhs[8]);
    Jii = mxGetScalar(prhs[9]);

    Kex = (int)mxGetScalar(prhs[10]);
    Kix = (int)mxGetScalar(prhs[11]);

    Kee = (int)mxGetScalar(prhs[12]);
    Kei = (int)mxGetScalar(prhs[13]);
    Kie = (int)mxGetScalar(prhs[14]);
    Kii = (int)mxGetScalar(prhs[15]);

    betaex = mxGetScalar(prhs[16]);
    betaix = mxGetScalar(prhs[17]);

    betaee = mxGetScalar(prhs[18]);
    betaei = mxGetScalar(prhs[19]);
    betaie = mxGetScalar(prhs[20]);
    betaii = mxGetScalar(prhs[21]);

    C = mxGetPr(prhs[22]);
    m1 = mxGetN(prhs[22]);
    m2 = mxGetM(prhs[22]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");
    gl = mxGetPr(prhs[23]);
    m1 = mxGetN(prhs[23]);
    m2 = mxGetM(prhs[23]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");
    Vleak = mxGetPr(prhs[24]);
    m1 = mxGetN(prhs[24]);
    m2 = mxGetM(prhs[24]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");
    DeltaT = mxGetPr(prhs[25]);
    m1 = mxGetN(prhs[25]);
    m2 = mxGetM(prhs[25]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");
    VT = mxGetPr(prhs[26]);
    m1 = mxGetN(prhs[26]);
    m2 = mxGetM(prhs[26]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");
    tref = mxGetPr(prhs[27]);
    m1 = mxGetN(prhs[27]);
    m2 = mxGetM(prhs[27]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");
    Vth = mxGetPr(prhs[28]);
    m1 = mxGetN(prhs[28]);
    m2 = mxGetM(prhs[28]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");
    Vre = mxGetPr(prhs[29]);
    m1 = mxGetN(prhs[29]);
    m2 = mxGetM(prhs[29]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");
    Vlb = mxGetPr(prhs[30]);
    m1 = mxGetN(prhs[30]);
    m2 = mxGetM(prhs[30]);
    if (m1 * m2 != 2)
        mexErrMsgTxt("All neuron parameters should be 2x1");

    tausynx = mxGetScalar(prhs[31]);
    tausyne = mxGetScalar(prhs[32]);
    tausyni = mxGetScalar(prhs[33]);

    v0 = mxGetPr(prhs[34]);
    N = mxGetM(prhs[34]);
    m2 = mxGetN(prhs[34]);
    if (N == 1 && m2 != 1)
        N = m2;

    T = mxGetScalar(prhs[35]);
    dt = mxGetScalar(prhs[36]);

    maxns = ((int)mxGetScalar(prhs[37]));

    Irecord = mxGetPr(prhs[38]);
    Nrecord = mxGetN(prhs[38]);

    m2 = mxGetM(prhs[38]);
    if (m2 != 2)
        mexErrMsgTxt("Irecord should be Nx2.");

    Irecord_I = mxGetPr(prhs[39]);
    Nrecord_I = mxGetN(prhs[39]);

    Irecord_I = mxGetPr(prhs[39]);
    Nrecord_I = mxGetN(prhs[39]);

    m2 = mxGetM(prhs[39]);
    if (m2 != 2)
        mexErrMsgTxt("Irecord should be Nx2.");

    /******
     * Finished importing variables.
     *******/

    /* Total number of each type of neuron */
    Ne = Ne1 * Ne1;
    Ni = Ni1 * Ni1;
    Nx = Nx1 * Nx1;

    /* Check for consistency with total number of neurons */
    if (N != Ne + Ni)
        mexErrMsgTxt("Ne1 and/or Ni1 not consistent with size of V0");

    /* Numebr of time bins */
    Nt = (int)(T / dt);

    /******
     * Now allocate output variables and temporary arrays.
     * This is also mex specific.  Use malloc in C, etc.
     *****/

    /* Allocate output vector */
    plhs[0] = mxCreateDoubleMatrix(3, maxns, mxREAL);
    s = mxGetPr(plhs[0]);

    /* Allocate output vector */
    plhs[1] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaxr = mxGetPr(plhs[1]);

    /* Allocate output vector */
    plhs[2] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaer = mxGetPr(plhs[2]);

    /* Allocate output vector */
    plhs[3] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphair = mxGetPr(plhs[3]);

    /* Allocate output vector */
    plhs[4] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    vr = mxGetPr(plhs[4]);

    plhs[5] = mxCreateDoubleMatrix(N, 1, mxREAL); // Initialize the scalar with a value, e.g., 0.0
    total_ffwd_input = mxGetPr(plhs[5]);          // Get a pointer to the scalar

    plhs[6] = mxCreateDoubleMatrix(N, 1, mxREAL); // Initialize the scalar for total_recurrent_ex_input
    total_recurrent_ex_input = mxGetPr(plhs[6]);

    plhs[7] = mxCreateDoubleMatrix(N, 1, mxREAL); // Initialize the scalar for total_recurrent_inh_input
    total_recurrent_inh_input = mxGetPr(plhs[7]);

    plhs[8] = mxCreateDoubleScalar(0.0); // Initialize the scalar for total_input
    total_input = mxGetPr(plhs[8]);

    plhs[9] = mxCreateDoubleScalar(0.0);
    input_count = mxGetPr(plhs[9]);

    plhs[10] = mxCreateDoubleMatrix(N, 1, mxREAL);
    input = mxGetPr(plhs[10]);

    plhs[11] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaer1 = mxGetPr(plhs[11]);

    plhs[12] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphair1 = mxGetPr(plhs[12]);

    plhs[13] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaxr1 = mxGetPr(plhs[13]);

    plhs[14] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaer2 = mxGetPr(plhs[14]);

    plhs[15] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphair2 = mxGetPr(plhs[15]);

    plhs[16] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaxr2 = mxGetPr(plhs[16]);

    plhs[17] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaer3 = mxGetPr(plhs[17]);

    plhs[18] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphair3 = mxGetPr(plhs[18]);

    plhs[19] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaxr3 = mxGetPr(plhs[19]);

    plhs[20] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaer4 = mxGetPr(plhs[20]);

    plhs[21] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphair4 = mxGetPr(plhs[21]);

    plhs[22] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    alphaxr4 = mxGetPr(plhs[22]);

    plhs[23] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    input_NL = mxGetPr(plhs[23]);

    plhs[24] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    input_NL_after = mxGetPr(plhs[24]);

    /* Allocate membrane potential */
    temp0 = mxCreateDoubleMatrix(N, 1, mxREAL);
    v = mxGetPr(temp0);

    temp1 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextE = mxGetPr(temp1);

    temp2 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextI = mxGetPr(temp2);

    temp3 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphae = mxGetPr(temp3);

    temp4 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphai = mxGetPr(temp4);

    /* Allocate membrane potential and other variables */

    // Allocate memory for JnextE1, JnextI1, and related variables
    temp5 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextE1 = mxGetPr(temp5);

    temp6 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextI1 = mxGetPr(temp6);

    temp7 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphae1 = mxGetPr(temp7);

    temp8 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphai1 = mxGetPr(temp8);

    // Allocate memory for JnextE2, JnextI2, and related variables
    temp9 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextE2 = mxGetPr(temp9);

    temp10 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextI2 = mxGetPr(temp10);

    temp11 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphae2 = mxGetPr(temp11);

    temp12 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphai2 = mxGetPr(temp12);

    // Allocate memory for JnextE3, JnextI3, and related variables
    temp13 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextE3 = mxGetPr(temp13);

    temp14 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextI3 = mxGetPr(temp14);

    temp15 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphae3 = mxGetPr(temp15);

    temp16 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphai3 = mxGetPr(temp16);

    // Allocate memory for JnextE4, JnextI4, and related variables
    temp17 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextE4 = mxGetPr(temp17);

    temp18 = mxCreateDoubleMatrix(N, 1, mxREAL);
    JnextI4 = mxGetPr(temp18);

    temp19 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphae4 = mxGetPr(temp19);

    temp20 = mxCreateDoubleMatrix(N, 1, mxREAL);
    alphai4 = mxGetPr(temp20);

    temp21 = mxCreateDoubleMatrix(N, 1, mxREAL);
    input_dendrite_1_before = mxGetPr(temp21);

    temp22 = mxCreateDoubleMatrix(N, 1, mxREAL);
    input_dendrite_1_after = mxGetPr(temp22);

    JnextX = mxMalloc(N * sizeof(double));
    JnextX1 = mxMalloc(N * sizeof(double));
    JnextX2 = mxMalloc(N * sizeof(double));
    JnextX3 = mxMalloc(N * sizeof(double));
    JnextX4 = mxMalloc(N * sizeof(double));
    alphax = mxMalloc(N * sizeof(double));
    alphax1 = mxMalloc(N * sizeof(double));
    alphax2 = mxMalloc(N * sizeof(double));
    alphax3 = mxMalloc(N * sizeof(double));
    alphax4 = mxMalloc(N * sizeof(double));

    /* Vectors for postsynaptic connections */
    Wee1 = mxMalloc(Ne * Kee * sizeof(int));
    Wee2 = mxMalloc(Ne * Kee * sizeof(int));
    Wei1 = mxMalloc(Ni * Kei * sizeof(int));
    Wei2 = mxMalloc(Ni * Kei * sizeof(int));
    Wie1 = mxMalloc(Ne * Kie * sizeof(int));
    Wie2 = mxMalloc(Ne * Kie * sizeof(int));
    Wii1 = mxMalloc(Ni * Kii * sizeof(int));
    Wii2 = mxMalloc(Ni * Kii * sizeof(int));
    Wex1 = mxMalloc(Nx * Kex * sizeof(int));
    Wex2 = mxMalloc(Nx * Kex * sizeof(int));
    Wix1 = mxMalloc(Nx * Kix * sizeof(int));
    Wix2 = mxMalloc(Nx * Kix * sizeof(int));
    Dee = mxMalloc(Ne * Kee * sizeof(int));
    Dei = mxMalloc(Ni * Kei * sizeof(int));
    Die = mxMalloc(Ne * Kie * sizeof(int));
    Dii = mxMalloc(Ni * Kii * sizeof(int));
    Dex = mxMalloc(Nx * Kex * sizeof(int));
    Dix = mxMalloc(Nx * Kix * sizeof(int));

    tempWee = mxMalloc(Kee * sizeof(int));

    tempWei = mxMalloc(Kei * sizeof(int));
    tempWie = mxMalloc(Kie * sizeof(int));
    tempWii = mxMalloc(Kii * sizeof(int));
    tempWex = mxMalloc(Kex * sizeof(int));
    tempWix = mxMalloc(Kix * sizeof(int));
    refstate = mxMalloc(N * sizeof(int));

    /*****
     * Finished allocating variables
     ****/

    int dendrite_type[N];
    int ffwd_dendrites[Nx];

    // double total_ffwd_input = 0.0;
    // double total_recurrent_ex_input = 0.0;
    // double total_recurrent_inh_input = 0.0;
    // *total_ffwd_input = 0.0;
    // *total_recurrent_ex_input = 0.0;
    // *total_recurrent_inh_input = 0.0;
    *total_input = 0.0;
    *input_count = 0.0;

    /* Inititalize variables */
    for (j = 0; j < N; j++)
    {
        v[j] = v0[j];
        refstate[j] = 0;
        JnextE[j] = 0;
        JnextI[j] = 0;
        alphae[j] = 0;
        alphai[j] = 0;
        alphax[j] = 0;
        JnextX[j] = 0;
        JnextX1[j] = 0;
        JnextX2[j] = 0;
        JnextX3[j] = 0;
        JnextX4[j] = 0;
        dendrite_type[j] = 0;
        input[j] = 0;
        total_ffwd_input[j] = 0;
        total_recurrent_ex_input[j] = 0;
        total_recurrent_inh_input[j] = 0;
        JnextE1[j] = 0;
        JnextI1[j] = 0;
        alphae1[j] = 0;
        alphai1[j] = 0;
        JnextE2[j] = 0;
        JnextI2[j] = 0;
        alphae2[j] = 0;
        alphai2[j] = 0;
        JnextE3[j] = 0;
        JnextI3[j] = 0;
        alphae3[j] = 0;
        alphai3[j] = 0;
        JnextE4[j] = 0;
        JnextI4[j] = 0;
        alphae4[j] = 0;
        alphai4[j] = 0;
        alphax1[j] = 0;
        alphax2[j] = 0;
        alphax3[j] = 0;
        alphax4[j] = 0;
        input_dendrite_1_before[j] = 0;
        input_dendrite_1_after[j] = 0;
    }

    mexPrintf("\nBuilding Network Architecure\n");
    mexEvalString("drawnow;");

    /* Record input currents and membrane potentials at first time bin */
    for (jj = 0; jj < Nrecord; jj++)
    {
        /* Find index into local variables */

        /* Find index into local variables */
        j1 = (int)round(Irecord[2 * jj] - 1);
        j2 = (int)round(Irecord[2 * jj + 1] - 1);

        /* Convert to 1D index, j */
        if (j1 < Ne1 && j2 < Ne1)
        {
            j = j1 + Ne1 * j2;
        }
        else if (j1 >= Ne1 && j2 >= Ne1)
        {
            j = (j1 - Ne1) + (j2 - Ne1) * Ni1 + Ne;
        }
        else
            mexErrMsgTxt("Indices in Irecord must have both terms <Ne1 or both terms >Ne1");

        if (j >= N || j < 0)
        {
            mexErrMsgTxt("Irecord contains out of bounds indices.");
        }

        /* Store currents and membrane potentials at first time bin */
        alphaer[jj + Nrecord * 0] = alphae[j];
        alphair[jj + Nrecord * 0] = alphai[j];
        alphaxr[jj + Nrecord * 0] = alphax[j];

        alphaer1[jj + Nrecord * 0] = alphae1[j];
        alphair1[jj + Nrecord * 0] = alphai1[j];
        alphaxr1[jj + Nrecord * 0] = alphax1[j];

        alphaer2[jj + Nrecord * 0] = alphae2[j];
        alphair2[jj + Nrecord * 0] = alphai2[j];
        alphaxr2[jj + Nrecord * 0] = alphax2[j];

        alphaer3[jj + Nrecord * 0] = alphae3[j];
        alphair3[jj + Nrecord * 0] = alphai3[j];
        alphaxr3[jj + Nrecord * 0] = alphax3[j];

        alphaer4[jj + Nrecord * 0] = alphae4[j];
        alphair4[jj + Nrecord * 0] = alphai4[j];
        alphaxr4[jj + Nrecord * 0] = alphax4[j];

        input_NL[jj + Nrecord * 0] = input_dendrite_1_before[j];
        input_NL_after[jj + Nrecord * 0] = input_dendrite_1_after[j];

        vr[jj + Nrecord * 0] = v[j];
    }

    /* Initialize connections */
    for (j = 0; j < Ne; j++)
    {

        /* Find index of cell along each dimension */
        j1 = j / Ne1;
        j2 = j % Ne1;

        /* Generate vectors of exc and inh postsynaptic targets (from excitatory neurons)*/
        CircRandNfun(tempWee, (double)j1, betaee, 0, Ne1 - 1, Kee);
        for (k = 0; k < Kee; k++)
            Wee1[j * Kee + k] = tempWee[k];
        CircRandNfun(tempWee, (double)j2, betaee, 0, Ne1 - 1, Kee);
        for (k = 0; k < Kee; k++)
        {
            Wee2[j * Kee + k] = tempWee[k];
            Dee[j * Kee + k] = closest_dendrite(j1, j2, Wee1[j * Kee + k], Wee2[j * Kee + k], Ne1 - 1, Ne1 - 1);
        }
        save1DArrayToCSV("ee_connections.csv", Dee, Kee);
        CircRandNfun(tempWie, ((double)j1 * ((double)Ni1 / (double)Ne1)), betaie, 0, Ni1 - 1, Kie);
        for (k = 0; k < Kie; k++)
            Wie1[j * Kie + k] = tempWie[k];
        CircRandNfun(tempWie, ((double)j2 * ((double)Ni1 / (double)Ne1)), betaie, 0, Ni1 - 1, Kie);
        for (k = 0; k < Kie; k++)
        {
            Wie2[j * Kie + k] = tempWie[k];
            Die[j * Kie + k] = closest_dendrite(((double)j1 * ((double)Ni1 / (double)Ne1)), ((double)j2 * ((double)Ni1 / (double)Ne1)), Wie1[j * Kie + k], Wie2[j * Kie + k], ((double)(Ne1 - 1) * ((double)Ni1 / (double)Ne1)), ((double)(Ne1 - 1) * ((double)Ni1 / (double)Ne1)));
        }
        save1DArrayToCSV("ie_connections.csv", Die, Kie);
    }

    for (j = 0; j < Ni; j++)
    {

        /* Find index of cell along each dimension */
        j1 = j / Ni1;
        j2 = j % Ni1;

        /* Generate vectors of exc and inh postsynaptic targets (from inhibitory neurons)*/
        CircRandNfun(tempWei, ((double)j1 * ((double)Ne1 / (double)Ni1)), betaei, 0, Ne1 - 1, Kei);
        for (k = 0; k < Kei; k++)
            Wei1[j * Kei + k] = tempWei[k];
        CircRandNfun(tempWei, ((double)j2 * ((double)Ne1 / (double)Ni1)), betaei, 0, Ne1 - 1, Kei);
        for (k = 0; k < Kei; k++)
        {
            Wei2[j * Kei + k] = tempWei[k];
            Dei[j * Kei + k] = closest_dendrite(((double)j1 * ((double)Ne1 / (double)Ni1)), ((double)j2 * ((double)Ne1 / (double)Ni1)), Wei1[j * Kei + k], Wei2[j * Kei + k], ((double)(Ni1 - 1) * ((double)Ne1 / (double)Ni1)), ((double)(Ni1 - 1) * ((double)Ne1 / (double)Ni1)));
        }
        save1DArrayToCSV("ei_connections.csv", Dei, Kei);
        CircRandNfun(tempWii, (double)j1, betaii, 0, Ni1 - 1, Kii);
        for (k = 0; k < Kii; k++)
            Wii1[j * Kii + k] = tempWii[k];
        CircRandNfun(tempWii, (double)j2, betaii, 0, Ni1 - 1, Kii);
        for (k = 0; k < Kii; k++)
        {
            Wii2[j * Kii + k] = tempWii[k];
            Dii[j * Kii + k] = closest_dendrite(j1, j2, Wii1[j * Kii + k], Wii2[j * Kii + k], Ni1 - 1, Ni1 - 1);
        }
        save1DArrayToCSV("ii_connections.csv", Dii, Kii);
    }

    for (j = 0; j < Nx; j++)
    {

        /* Find index of cell along each dimension */
        j1 = j / Nx1;
        j2 = j % Nx1;

        /* Generate vectors of exc and inh postsynaptic targets (from feedforward layer))*/
        CircRandNfun(tempWex, ((double)j1 * ((double)Ne1 / (double)Nx1)), betaex, 0, Ne1 - 1, Kex);
        for (k = 0; k < Kex; k++)
            Wex1[j * Kex + k] = tempWex[k];
        CircRandNfun(tempWex, ((double)j2 * ((double)Ne1 / (double)Nx1)), betaex, 0, Ne1 - 1, Kex);
        for (k = 0; k < Kex; k++)
        {
            Wex2[j * Kex + k] = tempWex[k];
            Dex[j * Kex + k] = closest_dendrite(((double)j1 * ((double)Ne1 / (double)Nx1)), ((double)j2 * ((double)Ne1 / (double)Nx1)), Wex1[j * Kex + k], Wex2[j * Kex + k], ((double)(Nx1 - 1) * ((double)Ne1 / (double)Nx1)), ((double)(Nx1 - 1) * ((double)Ne1 / (double)Nx1)));
        }
        save1DArrayToCSV("ex_connections.csv", Dex, Kex);
        CircRandNfun(tempWix, ((double)j1 * ((double)Ni1 / (double)Nx1)), betaix, 0, Ni1 - 1, Kix);
        for (k = 0; k < Kix; k++)
            Wix1[j * Kix + k] = tempWix[k];
        CircRandNfun(tempWix, ((double)j2 * ((double)Ni1 / (double)Nx1)), betaix, 0, Ni1 - 1, Kix);
        for (k = 0; k < Kix; k++)
        {
            Wix2[j * Kix + k] = tempWix[k];
            Dix[j * Kix + k] = closest_dendrite(((double)j1 * ((double)Ni1 / (double)Nx1)), ((double)j2 * ((double)Ni1 / (double)Nx1)), Wix1[j * Kix + k], Wix2[j * Kix + k], ((double)(Nx1 - 1) * ((double)Ni1 / (double)Nx1)), ((double)(Nx1 - 1) * ((double)Ni1 / (double)Nx1)));
        }
        save1DArrayToCSV("ix_connections.csv", Dix, Kix);
    }

    Ntref[0] = (int)round(tref[0] / dt);
    Ntref[1] = (int)round(tref[1] / dt);

    printfperiod = (int)(round(Nt / 20.0));

    /* Initialize number of spikes */
    ns = 0;

    iXspike = 0;

    /* Print portion complete every printperiod steps */
    mexPrintf("\n%d\n", printfperiod);
    mexEvalString("drawnow;");

    /* Main loop */
    /* Exit loop and issue a warning if max number of spikes is exceeded */
    for (i = 1; i < Nt && ns < maxns; i++)
    {

        /* Find all spikes in feedforward layer at this time bin */
        /* Add to corresponding elements of JnextX */
        while (sx[iXspike * 3 + 0] <= i * dt && iXspike < Nsx)
        {

            /* Find the index of each neuron in ffwd layer that spiked */
            int x_position = (int)round(sx[iXspike * 3 + 1] - 1); // X coordinate of presyn cell
            int y_position = (int)round(sx[iXspike * 3 + 2] - 1); // Y coordinate of presyn cell

            jspike = x_position * Nx1 + y_position;

            if (jspike < 0 || jspike >= Nx)
            {
                mexErrMsgTxt("Out of bounds index in sx.");
            }
            for (k = 0; k < Kex; k++)
            {
                postcell = Wex1[jspike * Kex + k] * Ne1 + Wex2[jspike * Kex + k];
                if (postcell < 0 || postcell >= N)
                    mexErrMsgTxt("Out of bounds index ino JnextX");
                JnextX[postcell] += Jex;
                int closest_dendrite = Dex[jspike * Kex + k];
                if (closest_dendrite == 0)
                {
                    JnextX1[postcell] += Jex;
                }
                else if (closest_dendrite == 1)
                {
                    JnextX2[postcell] += Jex;
                }
                else if (closest_dendrite == 2)
                {
                    JnextX3[postcell] += Jex;
                }
                else if (closest_dendrite == 3)
                {
                    JnextX4[postcell] += Jex;
                }
            }
            for (k = 0; k < Kix; k++)
            {
                postcell = Ne + Wix1[jspike * Kix + k] * Ni1 + Wix2[jspike * Kix + k];
                if (postcell < 0 || postcell >= N)
                    mexErrMsgTxt("Out of bounds index ino JnextX i");
                JnextX[postcell] += Jix;
                int closest_dendrite = Dix[jspike * Kix + k];
                if (closest_dendrite == 0)
                {
                    JnextX1[postcell] += Jix;
                }
                else if (closest_dendrite == 1)
                {
                    JnextX2[postcell] += Jix;
                }
                else if (closest_dendrite == 2)
                {
                    JnextX3[postcell] += Jix;
                }
                else if (closest_dendrite == 3)
                {
                    JnextX4[postcell] += Jix;
                }
            }

            iXspike++;
        }

        for (j = 0; j < N; j++)
        {
            /* Update synaptic variables */
            alphae[j] -= alphae[j] * (dt / tausyne);
            alphai[j] -= alphai[j] * (dt / tausyni);
            alphax[j] -= alphax[j] * (dt / tausynx);

            alphae1[j] -= alphae1[j] * (dt / tausyne);
            alphai1[j] -= alphai1[j] * (dt / tausyni);
            alphax1[j] -= alphax1[j] * (dt / tausynx);
            alphae2[j] -= alphae2[j] * (dt / tausyne);
            alphai2[j] -= alphai2[j] * (dt / tausyni);
            alphax2[j] -= alphax2[j] * (dt / tausynx);
            alphae3[j] -= alphae3[j] * (dt / tausyne);
            alphai3[j] -= alphai3[j] * (dt / tausyni);
            alphax3[j] -= alphax3[j] * (dt / tausynx);
            alphae4[j] -= alphae4[j] * (dt / tausyne);
            alphai4[j] -= alphai4[j] * (dt / tausyni);
            alphax4[j] -= alphax4[j] * (dt / tausynx);

            total_ffwd_input[j] = alphax[j];
            total_recurrent_ex_input[j] = alphae[j];
            total_recurrent_inh_input[j] = alphai[j];
            *total_input += alphax[j] + alphae[j] + alphai[j];

            input[j] = alphax[j] + alphae[j] + alphai[j];
            // mexPrintf("total input: %f\n", input);
            *input_count += 1.0;

            double total_input_current = alphae[j] + alphai[j] + alphax[j];
            double average_input_per_dendrite = 0.0;

            double noise = (double)rand() / RAND_MAX * 4.0f - 2.0f;
            noise = 0;

            double leak_current = (gl[0] * (v[j] - Vleak[0]));
            double leak_current_per_dendrite = leak_current / 4.0f;
            // leak_current = 0.0;

            double total_input_dendrite_1 = alphae1[j] + alphai1[j] + alphax1[j] + noise - leak_current_per_dendrite;
            double total_input_dendrite_2 = alphae2[j] + alphai2[j] + alphax2[j] + noise - leak_current_per_dendrite;
            double total_input_dendrite_3 = alphae3[j] + alphai3[j] + alphax3[j] + noise - leak_current_per_dendrite;
            double total_input_dendrite_4 = alphae4[j] + alphai4[j] + alphax4[j] + noise - leak_current_per_dendrite;

            input_dendrite_1_before[j] = total_input_dendrite_1;
            input_dendrite_1_after[j] = total_input_dendrite_1;

            double centered_sigmoid2 = total_input_dendrite_1 + total_input_dendrite_2 + total_input_dendrite_3 + total_input_dendrite_4;

            if (j < Ne)
            { /* If cell is excitatory */

                if (refstate[j] <= 0)
                    v[j] += fmax(centered_sigmoid2 * dt / C[0], Vlb[0] - v[j]);
                else
                {
                    if (refstate[j] > 1)
                        v[j] = Vth[0];
                    else
                        v[j] = Vre[0];
                    refstate[j]--;
                }

                /* If a spike occurs */
                if (v[j] >= Vth[0] && refstate[j] <= 0 && ns < maxns)
                {
                    refstate[j] = Ntref[0];
                    v[j] = Vth[0];               /* reset membrane potential */
                    s[0 + 3 * ns] = i * dt;      /* spike time */
                    s[2 + 3 * ns] = j / Ne1 + 1; /* neuron index 1 (row?) */
                    s[1 + 3 * ns] = j % Ne1 + 1; /* neuron index 2 */

                    ns++; /* update total number of spikes */

                    // int pre_x = j / Ne1; /* Pre-synaptic row (y-coordinate) */
                    // int pre_y = j % Ne1; /* Pre-synaptic column (x-coordinate) */

                    /* For each postsynaptic target, propagate spike into JnextE */
                    for (k = 0; k < Kee; k++)
                    {
                        int postsynaptic_cell = Wee1[j * Kee + k] * Ne1 + Wee2[j * Kee + k];
                        JnextE[postsynaptic_cell] += Jee;
                        int closest_dendrite = Dee[j * Kee + k];
                        if (closest_dendrite == 0)
                        {
                            JnextE1[postsynaptic_cell] += Jee;
                        }
                        else if (closest_dendrite == 1)
                        {
                            JnextE2[postsynaptic_cell] += Jee;
                        }
                        else if (closest_dendrite == 2)
                        {
                            JnextE3[postsynaptic_cell] += Jee;
                        }
                        else if (closest_dendrite == 3)
                        {
                            JnextE4[postsynaptic_cell] += Jee;
                        }
                    }
                    for (k = 0; k < Kie; k++)
                    {
                        int postsynaptic_cell = Ne + Wie1[j * Kie + k] * Ni1 + Wie2[j * Kie + k];
                        JnextE[postsynaptic_cell] += Jie;
                        int closest_dendrite = Die[j * Kie + k];
                        if (closest_dendrite == 0)
                        {
                            JnextE1[postsynaptic_cell] += Jie;
                        }
                        else if (closest_dendrite == 1)
                        {
                            JnextE2[postsynaptic_cell] += Jie;
                        }
                        else if (closest_dendrite == 2)
                        {
                            JnextE3[postsynaptic_cell] += Jie;
                        }
                        else if (closest_dendrite == 3)
                        {
                            JnextE4[postsynaptic_cell] += Jie;
                        }
                    }
                }
            }

            else
            { /* If cell is inhibitory */

                if (refstate[j] <= 0)
                    v[j] += fmax(centered_sigmoid2 * dt / C[1], Vlb[1] - v[j]);
                else
                {
                    if (refstate[j] > 1)
                        v[j] = Vth[1];
                    else
                        v[j] = Vre[1];
                    refstate[j]--;
                }

                /* If a spike occurs */
                if (v[j] >= Vth[1] && refstate[j] <= 0 && ns < maxns)
                {
                    refstate[j] = Ntref[1];
                    v[j] = Vth[1];                         /* reset membrane potential */
                    s[0 + 3 * ns] = i * dt;                /* spike time */
                    s[2 + 3 * ns] = -((j - Ne) / Ni1) - 1; /* neuron index 1 */
                    s[1 + 3 * ns] = -((j - Ne) % Ni1) - 1; /* neuron index 2 */

                    ns++; /* update total number of spikes */

                    // int pre_x = (j - Ne) / Ni1; /* Pre-synaptic row (y-coordinate) */
                    // int pre_y = (j - Ne) % Ni1; /* Pre-synaptic column (x-coordinate) */

                    /* For each postsynaptic target, propagate spike into JnextI */
                    for (k = 0; k < Kei; k++)
                    {
                        int postsynaptic_cell = Wei1[(j - Ne) * Kei + k] * Ne1 + Wei2[(j - Ne) * Kei + k];
                        JnextI[postsynaptic_cell] += Jei;
                        int closest_dendrite = Dei[(j - Ne) * Kei + k];
                        if (closest_dendrite == 0)
                        {
                            JnextI1[postsynaptic_cell] += Jei;
                        }
                        else if (closest_dendrite == 1)
                        {
                            JnextI2[postsynaptic_cell] += Jei;
                        }
                        else if (closest_dendrite == 2)
                        {
                            JnextI3[postsynaptic_cell] += Jei;
                        }
                        else if (closest_dendrite == 3)
                        {
                            JnextI4[postsynaptic_cell] += Jei;
                        }
                    }
                    for (k = 0; k < Kii; k++)
                    {
                        int postsynaptic_cell = Ne + Wii1[(j - Ne) * Kii + k] * Ni1 + Wii2[(j - Ne) * Kii + k];
                        JnextI[postsynaptic_cell] += Jii;
                        int closest_dendrite = Dii[(j - Ne) * Kii + k];
                        if (closest_dendrite == 0)
                        {
                            JnextI1[postsynaptic_cell] += Jii;
                        }
                        else if (closest_dendrite == 1)
                        {
                            JnextI2[postsynaptic_cell] += Jii;
                        }
                        else if (closest_dendrite == 2)
                        {
                            JnextI3[postsynaptic_cell] += Jii;
                        }
                        else if (closest_dendrite == 3)
                        {
                            JnextI4[postsynaptic_cell] += Jii;
                        }
                    }
                }
            }
        }

        /* Store recorded variables */
        for (jj = 0; jj < Nrecord; jj++)
        {

            /* Find index into local variables */
            j1 = (int)round(Irecord[2 * jj + 0] - 1);
            j2 = (int)round(Irecord[2 * jj + 1] - 1);

            if (j1 < Ne1 && j2 < Ne1)
            {
                j = j1 + Ne1 * j2;
            }
            else if (j1 >= Ne1 && j2 >= Ne1)
            {
                j = (j1 - Ne1) + (j2 - Ne1) * Ni1 + Ne;
            }
            else
                mexErrMsgTxt("Indices in Irecord must have both terms <Ne1 or both terms >Ne1");

            alphaer[jj + Nrecord * i] = alphae[j];
            alphair[jj + Nrecord * i] = alphai[j];
            alphaxr[jj + Nrecord * i] = alphax[j];

            alphaer1[jj + Nrecord * i] = alphae1[j];
            alphair1[jj + Nrecord * i] = alphai1[j];
            alphaer2[jj + Nrecord * i] = alphae2[j];
            alphair2[jj + Nrecord * i] = alphai2[j];
            alphaer3[jj + Nrecord * i] = alphae3[j];
            alphair3[jj + Nrecord * i] = alphai3[j];
            alphaer4[jj + Nrecord * i] = alphae4[j];
            alphair4[jj + Nrecord * i] = alphai4[j];
            alphaxr1[jj + Nrecord * i] = alphax1[j];
            alphaxr2[jj + Nrecord * i] = alphax2[j];
            alphaxr3[jj + Nrecord * i] = alphax3[j];
            alphaxr4[jj + Nrecord * i] = alphax4[j];

            vr[jj + Nrecord * i] = v[j];

            input_NL[jj + Nrecord * i] = input_dendrite_1_before[j];
            input_NL_after[jj + Nrecord * i] = input_dendrite_1_after[j];
        }

        /* Use Jnext vectors to update synaptic variables */
        for (j = 0; j < N; j++)
        {
            // synaptic decay
            alphae[j] += JnextE[j] / tausyne;
            alphai[j] += JnextI[j] / tausyni;
            alphax[j] += JnextX[j] / tausynx;
            alphae1[j] += JnextE1[j] / tausyne;
            alphai1[j] += JnextI1[j] / tausyni;
            alphae2[j] += JnextE2[j] / tausyne;
            alphai2[j] += JnextI2[j] / tausyni;
            alphae3[j] += JnextE3[j] / tausyne;
            alphai3[j] += JnextI3[j] / tausyni;
            alphae4[j] += JnextE4[j] / tausyne;
            alphai4[j] += JnextI4[j] / tausyni;
            alphax1[j] += JnextX1[j] / tausynx;
            alphax2[j] += JnextX2[j] / tausynx;
            alphax3[j] += JnextX3[j] / tausynx;
            alphax4[j] += JnextX4[j] / tausynx;
            JnextE[j] = 0;
            JnextI[j] = 0;
            JnextX[j] = 0;
            JnextE1[j] = 0;
            JnextI1[j] = 0;
            JnextE2[j] = 0;
            JnextI2[j] = 0;
            JnextE3[j] = 0;
            JnextI3[j] = 0;
            JnextE4[j] = 0;
            JnextI4[j] = 0;
            JnextX1[j] = 0;
            JnextX2[j] = 0;
            JnextX3[j] = 0;
            JnextX4[j] = 0;
        }

        /* Print percent complete every printfperiod time steps
         * This might not actually print until the full simulation
         * is complete due to how some versions of Matlab treat the
         * drawnow signal coming from a mex file */
        if (i % printfperiod == 0)
        {
            mexPrintf("\n%d percent complete  rate = %2.2fHz", i * 100 / Nt, 1000 * ((double)(ns)) / (((double)(N)) * ((double)(i)) * dt));
            mexEvalString("drawnow;");
        }
    }

    /* Issue a warning if max number of spikes reached */
    if (ns >= maxns)
        mexWarnMsgTxt("Maximum number of spikes reached, simulation terminated.");

    mexPrintf("\nTotal ffwd input: %f\n", *total_ffwd_input);
    mexPrintf("\nTotal recurrent exc input: %f\n", *total_recurrent_ex_input);
    mexPrintf("\nTotal recurrent inh input: %f\n", *total_recurrent_inh_input);

    /* Free allocated memory */
    mxDestroyArray(temp0);
    mxDestroyArray(temp1);
    mxDestroyArray(temp2);
    mxDestroyArray(temp3);
    mxDestroyArray(temp4);
    mxDestroyArray(temp5);
    mxDestroyArray(temp6);
    mxDestroyArray(temp7);
    mxDestroyArray(temp8);
    mxDestroyArray(temp9);
    mxDestroyArray(temp10);
    mxDestroyArray(temp11);
    mxDestroyArray(temp12);
    mxDestroyArray(temp13);
    mxDestroyArray(temp14);
    mxDestroyArray(temp15);
    mxDestroyArray(temp16);
    mxDestroyArray(temp17);
    mxDestroyArray(temp18);
    mxDestroyArray(temp19);
    mxDestroyArray(temp20);
    mxFree(Wee1);
    mxFree(Wee2);
    mxFree(Wei1);
    mxFree(Wei2);
    mxFree(Wie1);
    mxFree(Wie2);
    mxFree(Wii1);
    mxFree(Wii2);
    mxFree(Wex1);
    mxFree(Wex2);
    mxFree(Wix1);
    mxFree(Wix2);
    mxFree(tempWee);
    mxFree(tempWei);
    mxFree(tempWie);
    mxFree(tempWii);
    mxFree(tempWix);
    mxFree(tempWex);
    mxFree(JnextX);
    mxFree(JnextX1);
    mxFree(JnextX2);
    mxFree(JnextX3);
    mxFree(JnextX4);
    mxFree(alphax);
    mxFree(alphax1);
    mxFree(alphax2);
    mxFree(alphax3);
    mxFree(alphax4);
    mxFree(refstate);
    mxFree(Dex);
    mxFree(Dix);
    mxFree(Dee);
    mxFree(Dei);
    mxFree(Die);
    mxFree(Dii);
}

/* This function generates n random variables that are integers between min and max.
 * The distribution of these random variables is like a Gaussian distribution with
 * mean mu and std sigma, but rounded to the nearest integer and wrapped around the
 * interval [min,max]. The values are stored in the location pointed to by z, so you
 * better make sure that you allocated room for at least n integers in z.
 *
 *  creates random numbers that mimic a gaussian distribution with given meand and std, wrapped within defined interval [min, max]
 */
void CircRandNfun(int *z, double mu, double sigma, int min, int max, int n)
{
    int i;
    double u1, u2;
    double testy;
    int matlabmod(int, int);

    for (i = 1; i < n; i += 2)
    {
        /*
         * Generate two random numbers between 0 and 1 - uniformly distributed!!
         */
        u1 = drand48();
        u2 = drand48();

        /*
         * 1) Box-Muller transform to take uniformly distributed numbers to 2 independent normally distributed numbers
         *
         * z[i-1]=matlabmod((int)round(sigma*sqrt(-2*log(u1))*cos(2*M_PI*u2)+mu)-min,max-min+1)+min;
         * z[i]=matlabmod((int)round(sigma*sqrt(-2*log(u1))*sin(2*M_PI*u2)+mu)-min,max-min+1)+min;
         *
         * Generate a Gaussian random number using the Box-Muller transform
         *
         */
        double gaussian1 = sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2) + mu;
        double gaussian2 = sigma * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2) + mu;
        /*
         * Map the Gaussian result to the interval [min, max] using mod -> get distance in each direction from presyn neuron periodically
         */
        z[i - 1] = matlabmod((int)round(gaussian1) - min, max - min + 1) + min;
        z[i] = matlabmod((int)round(gaussian2) - min, max - min + 1) + min;
    }

    if (i == n)
    {
        u1 = drand48();
        u2 = drand48();
        z[i - 1] = (int)round(sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2) + mu);
        z[i - 1] = matlabmod(z[i - 1] - min, max - min + 1) + min;
    }
}

/* Implements a "mod" function that behaves like matlab's version
 * of mod instead of like % in C.  I forgot what the difference is,
 * but remember that it's important, especially when a is negative.
 */
int matlabmod(int a, int b)
{
    int c;
    c = a % b;
    while (c < 0)
        c += b;
    return c;
}

/**
 * Determines the closest dendrite to a given presynaptic location,
 *        accounting for periodic boundary conditions.
 *
 * This function computes the Euclidean distance from a presynaptic location
 * (pre_x, pre_y) to the four dendritic branches of a postsynaptic neuron
 * located at (post_x, post_y). It considers periodic boundary conditions in
 * both x and y directions. The function returns the index of the closest
 * dendrite. In the case of ties, it randomly selects one of the closest dendrites.
 *
 * pre_x The x-coordinate of the presynaptic neuron.
 * pre_y The y-coordinate of the presynaptic neuron.
 * post_x The x-coordinate of the postsynaptic neuron.
 * post_y The y-coordinate of the postsynaptic neuron.
 * Lx The width of the periodic domain in the x direction.
 * Ly The height of the periodic domain in the y direction.
 *
 * returns an integer representing the index of the closest dendrite:
 *         0 - Right dendrite
 *         1 - Left dendrite
 *         2 - Up dendrite
 *         3 - Down dendrite
 */
int closest_dendrite(double pre_x, double pre_y, double post_x, double post_y, int Lx, int Ly)
{

    double periodic_distance(double, double, int);
    double dendrite_length = 1.0;

    // Compute positions of each dendrite
    double dendrite_x_right = post_x + dendrite_length;
    if (dendrite_x_right > Lx)
    {
        dendrite_x_right = dendrite_x_right - Lx;
    }
    double dendrite_y_right = post_y;

    double dendrite_x_left = post_x - dendrite_length;
    if (dendrite_x_left < 0)
    {
        dendrite_x_left = dendrite_x_left + Lx;
    }
    double dendrite_y_left = post_y;

    double dendrite_x_up = post_x;
    double dendrite_y_up = post_y + dendrite_length;
    if (dendrite_y_up > Ly)
    {
        dendrite_y_up = dendrite_y_up - Ly;
    }

    double dendrite_x_down = post_x;
    double dendrite_y_down = post_y - dendrite_length;
    if (dendrite_y_down < 0)
    {
        dendrite_y_down = dendrite_y_down + Ly;
    }

    int dendrites = 4;
    double distances[dendrites];

    // Calculate distances accounting for periodic boundary conditions
    distances[0] = sqrt(pow(periodic_distance(pre_x, dendrite_x_right, Lx), 2) + pow(periodic_distance(pre_y, dendrite_y_right, Ly), 2));
    distances[1] = sqrt(pow(periodic_distance(pre_x, dendrite_x_left, Lx), 2) + pow(periodic_distance(pre_y, dendrite_y_left, Ly), 2));
    distances[2] = sqrt(pow(periodic_distance(pre_x, dendrite_x_up, Lx), 2) + pow(periodic_distance(pre_y, dendrite_y_up, Ly), 2));
    distances[3] = sqrt(pow(periodic_distance(pre_x, dendrite_x_down, Lx), 2) + pow(periodic_distance(pre_y, dendrite_y_down, Ly), 2));

    // Find the minimum distance, if some dendrites have the same distance, choose randomly
    int min_index = 0;
    double min_distance = distances[0];
    for (int i = 1; i < dendrites; i++)
    {
        if (distances[i] < min_distance)
        {
            min_distance = distances[i];
            min_index = i;
        }
        else if (distances[i] == min_distance)
        {
            if (rand() % 2 == 0)
            {
                min_distance = distances[i];
                min_index = i;
            }
        }
    }
    // print in 0.1 percent of teh cases:
    // if (rand() % 1000 == 0)
    // {
    //     mexPrintf("pre_x: %f, pre_y: %f, post_x: %f, post_y: %f, Lx: %d, Ly: %d\n", pre_x, pre_y, post_x, post_y, Lx, Ly);
    //     mexPrintf("dendrite_x_right: %f, dendrite_y_right: %f\n", dendrite_x_right, dendrite_y_right);
    //     mexPrintf("dendrite_x_left: %f, dendrite_y_left: %f\n", dendrite_x_left, dendrite_y_left);
    //     mexPrintf("dendrite_x_up: %f, dendrite_y_up: %f\n", dendrite_x_up, dendrite_y_up);
    //     mexPrintf("dendrite_x_down: %f, dendrite_y_down: %f\n", dendrite_x_down, dendrite_y_down);
    //     mexPrintf("distances[0]: %f, distances[1]: %f, distances[2]: %f, distances[3]: %f\n", distances[0], distances[1], distances[2], distances[3]);
    //     mexPrintf("min_distance: %f, min_index: %d\n", min_distance, min_index);
    // }

    return min_index;
}

// Function to compute the minimum periodic distance between two coordinates
// double periodic_distance(double coord1, double coord2, int L)
// {
//     double diff = fabs(coord1 - coord2);
//     double min = fmin(diff, L - diff);
//     return min; // Take the minimum of direct distance or wrapping around
// }
// def periodic_b_dendrite(x, length):
//     if x < 0:
//         x += length
//     elif x >= length:
//         x -= length
//     return x
double periodic_distance(double coord1, double coord2, int L)
{
    //   diff = abs(coord1 - coord2) % L
    //     return min(diff, L - diff)

    // return min(diff, L - diff)
    double diff = matlabmod(fabs(coord1 - coord2), L);
    double min = fmin(diff, L - diff);
    return min; // Take the minimum of direct distance or wrapping around
}

double applySigmoid(double input, double lambda, double a)
{

    // double a = 0.0;      // Center point of sigmoid

    // Apply sigmoid
    double sigmoid_input = 1.0 / (1.0 + exp(-lambda * (input - a))); // s(x)

    return sigmoid_input;
}

void save1DArrayToCSV(const char *filename, int *array, int length)
{
    // Open the file in write mode. This will create the file if it does not exist.
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < length; i++)
    {
        fprintf(file, "%d\n", array[i]);
    }

    fclose(file);
}