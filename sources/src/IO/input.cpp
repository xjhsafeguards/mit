#include "input.h"

auto INPUT = new Input;

void Input::init(){

}

void Input::iarg(int argc,char** argv){
  if(argc != 1)
  {
      for(int i=1; i<argc; ++i)
      {
          //if(strncmp(argv[i],"count",5) == 0 or strncmp(argv[i],"-c",2) == 0)//{INPUT.calculation="count";}
          //else if(strncmp(argv[i],"help",4) == 0 or strncmp(argv[i],"-h",2) == 0)
          //{Help::Routine(INPUT.calculation); return 1;}
          //else if(strncmp(argv[i],"--cal=",6) == 0)
          //{string tmp = argv[i];INPUT.calculation=tmp.substr(6);}
          //else if(strncmp(argv[i],"--type=",7) == 0)
          //{string tmp = argv[i];convstring(tmp.substr(7),INPUT.type);}
          //else if(strncmp(argv[i],"-t",2) == 0)
          //{i++; string tmp = argv[i];convstring(tmp,INPUT.type);}
          //else if(strncmp(argv[i],"--delta=",8) == 0)
          //{string tmp = argv[i];convstring(tmp.substr(8),INPUT.delta);}
          //else if(strncmp(argv[i],"-d",2) == 0)
          //{string tmp = argv[++i];convstring(tmp,INPUT.type);}
      }
  }
}
