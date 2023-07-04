function SturctToJSON(Cell,NameFile)
  Txt_JSON = jsonencode(Cell,'PrettyPrint',true) ;
  oid = fopen(NameFile,'wt');
  fprintf(oid, Txt_JSON);
  fclose(oid);
end
