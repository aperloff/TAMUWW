void list_files(const char *dirname="/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/TTH_HToBB_M125_JESUp/crab_run_MatrixElement_TTH_HToBB_M125_JESUp/150512_043119/0001/", const char *ext=".root")
{
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   int count = 0;
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(ext)) {
            cout << fname.Data() << endl;
            count++;
         }
      }
   }
   cout << "Total number of files: " << count << endl;
}
