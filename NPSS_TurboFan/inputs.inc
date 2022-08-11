// ---------------
//Ambient
// ---------------
Eng.Amb.alt_in            = 35000;
Eng.Amb.MN_in			  = 0.8;
// ---------------
//Fan
// ---------------
Eng.CmpF.PRdes            = 1.685;
Eng.CmpF.effDes           = 0.8948;
Eng.CmpF.Fl_O.MNdes       = 0.4578;
// ---------------
//Splitter 
// ---------------           
Eng.SpltFan.BPRdes        = 5.1205;
// ---------------
//LPC
// ---------------                    
Eng.CmpL.PRdes            = 1.935;
Eng.CmpL.effDes           = 0.9243;
Eng.CmpL.Fl_O.MNdes       = 0.3059;
// ---------------
// Interstage bleed
// ---------------
Eng.Core.CmpH.Fl_B_CH.fracBldW    = 0.050708;
Eng.Core.CmpH.Fl_B_CH.fracBldP    = 0.5;
Eng.Core.CmpH.Fl_B_CH.fracBldWork = 0.5;

Eng.Core.CmpH.Fl_B_NC.fracBldW    = 0.020274;
Eng.Core.CmpH.Fl_B_NC.fracBldP    = 0.55;
Eng.Core.CmpH.Fl_B_NC.fracBldWork = 0.5;

// Eng.Core.CmpH.Fl_B027.fracBldW    = 0.0445; // set in Core model based on des or offdes
Eng.Core.CmpH.Fl_B027.fracBldP    = 0.5;
Eng.Core.CmpH.Fl_B027.fracBldWork = 0.5;
// ---------------
//HPC
// ---------------
Eng.Core.CmpH.PRdes       = 9.369;
Eng.Core.CmpH.effDes      = 0.8723;
Eng.Core.CmpH.Fl_O.MNdes  = 0.2442;
// ---------------
//HPT Cooling 
// ---------------
Eng.Core.B030.TCLA_NC.fracW  = 0.067214;
Eng.Core.B030.TCLA_CH.fracW  = 0.101256;
// ---------------
//Burner
// ---------------
Eng.Core.BrnPri.dPqP_dmd  = 0.054;
// ---------------
//HPT
// ---------------
Eng.Core.TrbH.effDes      = 0.8888;
// ---------------
//LPT
// ---------------
Eng.TrbL.effDes           = 0.90;
// ---------------
//Core Nozzle
// ---------------
Eng.NozPri.Cv_in          = 0.9933;
// ---------------
//Byp Nozzle
// ---------------
Eng.NozSec.Cv_in          = 0.9939;


// RESET autosolversetup
Eng.Core.CmpH.ind_CmpH_eff.autoSetup = FALSE;
Eng.Core.dep_CmpH_polyEff.autoSetup  = FALSE;

Eng.CmpF.ind_CmpF_eff.autoSetup      = FALSE;
Eng.dep_CmpF_polyEff.autoSetup       = FALSE;

Eng.CmpL.ind_CmpL_eff.autoSetup      = FALSE;
Eng.dep_CmpL_polyEff.autoSetup       = FALSE;


