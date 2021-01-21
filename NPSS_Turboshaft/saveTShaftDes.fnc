OutFileStream file {
    filename = "Tshaft_design.int"; 
}

file<<"// LHV based on fuel\n";
file<<"Eng.FusEng.LHV   = "<<Eng.FusEng.LHV  <<";"<<endl;
file<<"Eng.FsEng.W_in   = "<<Eng.FsEng.W_in  <<";"<<endl;
file<<"Eng.FusEng.Wfuel = "<<Eng.FusEng.Wfuel<<";"<<endl<<endl;

file<<"// Map Scalars"<<endl;
file<<"Eng.CmpL.S_map.s_effDes = "<< Eng.CmpL.S_map.s_effDes <<";"<< endl;
file<<"Eng.CmpL.S_map.s_PRdes  = "<< Eng.CmpL.S_map.s_PRdes  <<";"<< endl;
file<<"Eng.CmpL.S_map.s_WcDes  = "<< Eng.CmpL.S_map.s_WcDes  <<";"<< endl;
file<<"Eng.CmpL.S_map.s_NcDes  = "<< Eng.CmpL.S_map.s_NcDes  <<";"<< endl<<endl;

file<<"Eng.CmpH.S_map.s_effDes = "<< Eng.CmpH.S_map.s_effDes <<";"<< endl;
file<<"Eng.CmpH.S_map.s_PRdes  = "<< Eng.CmpH.S_map.s_PRdes  <<";"<< endl;
file<<"Eng.CmpH.S_map.s_WcDes  = "<< Eng.CmpH.S_map.s_WcDes  <<";"<< endl;
file<<"Eng.CmpH.S_map.s_NcDes  = "<< Eng.CmpH.S_map.s_NcDes  <<";"<< endl<<endl;

file<<"Eng.TrbH.S_map.s_effDes = "<< Eng.TrbH.S_map.s_effDes <<";"<< endl;
file<<"Eng.TrbH.S_map.s_PRdes  = "<< Eng.TrbH.S_map.s_PRdes  <<";"<< endl;
file<<"Eng.TrbH.S_map.s_WpDes  = "<< Eng.TrbH.S_map.s_WpDes  <<";"<< endl;
file<<"Eng.TrbH.S_map.s_NpDes  = "<< Eng.TrbH.S_map.s_NpDes  <<";"<< endl<<endl;

file<<"Eng.TrbP.S_map.s_effDes = "<< Eng.TrbP.S_map.s_effDes <<";"<< endl;
file<<"Eng.TrbP.S_map.s_PRdes  = "<< Eng.TrbP.S_map.s_PRdes  <<";"<< endl;
file<<"Eng.TrbP.S_map.s_WpDes  = "<< Eng.TrbP.S_map.s_WpDes  <<";"<< endl;
file<<"Eng.TrbP.S_map.s_NpDes  = "<< Eng.TrbP.S_map.s_NpDes  <<";"<< endl<<endl;

file<< "Eng.CmpL.NcDes = "<< Eng.CmpL.NcDes<<";"<<endl;
file<< "Eng.CmpH.NcDes = "<< Eng.CmpH.NcDes<<";"<<endl;
file<< "Eng.TrbH.NpDes = "<< Eng.TrbH.NpDes<<";"<<endl;
file<< "Eng.TrbP.NpDes = "<< Eng.TrbP.NpDes<<";"<<endl<<endl;

file<< "Eng.CmpL.S_map.NcMapDes    = "<<Eng.CmpL.S_map.NcMapDes   <<";"<<endl;
file<< "Eng.CmpL.S_map.WcMapDes    = "<<Eng.CmpL.S_map.WcMapDes   <<";"<<endl;
file<< "Eng.CmpL.S_map.PRmapDes    = "<<Eng.CmpL.S_map.PRmapDes   <<";"<<endl;
file<< "Eng.CmpL.S_map.effMapDes   = "<<Eng.CmpL.S_map.effMapDes  <<";"<<endl;
file<< "Eng.CmpL.S_map.RlineMapDes = "<<Eng.CmpL.S_map.RlineMapDes<<";"<<endl;

file<< "Eng.CmpH.S_map.NcMapDes    = "<<Eng.CmpH.S_map.NcMapDes   <<";"<<endl;
file<< "Eng.CmpH.S_map.WcMapDes    = "<<Eng.CmpH.S_map.WcMapDes   <<";"<<endl;
file<< "Eng.CmpH.S_map.PRmapDes    = "<<Eng.CmpH.S_map.PRmapDes   <<";"<<endl;
file<< "Eng.CmpH.S_map.effMapDes   = "<<Eng.CmpH.S_map.effMapDes  <<";"<<endl;
file<< "Eng.CmpH.S_map.RlineMapDes = "<<Eng.CmpH.S_map.RlineMapDes<<";"<<endl;

file<< "Eng.TrbH.S_map.PRmapDes = "<<Eng.TrbH.S_map.PRmapDes <<";"<<endl;
file<< "Eng.TrbH.S_map.NpMapDes = "<<Eng.TrbH.S_map.NpMapDes <<";"<<endl;
file<< "Eng.TrbP.S_map.PRmapDes = "<<Eng.TrbP.S_map.PRmapDes <<";"<<endl;
file<< "Eng.TrbP.S_map.NpMapDes = "<<Eng.TrbP.S_map.NpMapDes <<";"<<endl;

//Nozzle
file<< "Eng.NozPri.AthCold = "<< Eng.NozPri.AthCold <<";"<< endl; 

//PCEC
file<<"Eng.PCEC.l    = "<<Eng.PCEC.l    <<";"<<endl;
file<<"Eng.PCEC.w    = "<<Eng.PCEC.w    <<";"<<endl;
file<<"Eng.PCEC.cpsi = "<<Eng.PCEC.cpsi <<";"<<endl;
file<<"Eng.PCEC.Af   = "<<Eng.PCEC.Af   <<";"<<endl<<endl;

file<<"Eng.B030.TCLA_NC.fracW = "<<Eng.B030.TCLA_NC.fracW <<";"<<endl;
file<<"Eng.B030.TCLA_CH.fracW = "<<Eng.B030.TCLA_CH.fracW <<";"<<endl<<endl;


file << "Eng.FsEng.Fl_O.Aphy  = "<<Eng.FsEng.Fl_O.Aphy  <<";"<<endl;
file << "Eng.InEng.Fl_O.Aphy  = "<<Eng.InEng.Fl_O.Aphy  <<";"<<endl;
file << "Eng.CmpL.Fl_O.Aphy   = "<<Eng.CmpL.Fl_O.Aphy   <<";"<<endl;
file << "Eng.D025.Fl_O.Aphy   = "<<Eng.D025.Fl_O.Aphy   <<";"<<endl;
file << "Eng.CmpH.Fl_O.Aphy   = "<<Eng.CmpH.Fl_O.Aphy   <<";"<<endl;
file << "Eng.B030.Fl_O.Aphy   = "<<Eng.B030.Fl_O.Aphy   <<";"<<endl;
file << "Eng.BrnPri.Fl_O.Aphy = "<<Eng.BrnPri.Fl_O.Aphy <<";"<<endl;
file << "Eng.TrbH.Fl_O.Aphy   = "<<Eng.TrbH.Fl_O.Aphy   <<";"<<endl;
file << "Eng.D050.Fl_O.Aphy   = "<<Eng.D050.Fl_O.Aphy   <<";"<<endl;
file << "Eng.TrbP.Fl_O.Aphy   = "<<Eng.TrbP.Fl_O.Aphy   <<";"<<endl;
file << "Eng.D060.Fl_O.Aphy   = "<<Eng.D060.Fl_O.Aphy   <<";"<<endl;
file << "Eng.PCEC.Fl_O.Aphy   = "<<Eng.PCEC.Fl_O.Aphy   <<";"<<endl;
file << "Eng.NozPri.Fl_O.Aphy = "<<Eng.NozPri.Fl_O.Aphy <<";"<<endl;