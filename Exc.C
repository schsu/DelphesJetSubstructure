void Exc () {
  gSystem->Load("JetSub");
  
  
  JetSub aaa ("./files/600.root", "files/test.root");
  aaa.Initialize();
  
  TStopwatch s;
  aaa.Run();
  cout << s.CpuTime() << endl;
  
  return;
}
