void failed(char *string){
     fprintf(stderr,"%s",string);
     exit(1);
}

/*********************************/

char *c1t (int n1)
{
  char *p;
  if ((p = (char *) malloc ((size_t) n1 * sizeof (char))) == NULL)
      failed ( (char *) "c1t: failed");
  return p;
}

/*********************************/

int **i2t (int n1, int n2)
{
  int **p, *a;
  int i;
  if ((p = (int **) malloc ((size_t) n1 * sizeof (int *))) == NULL)
      failed ( (char *) "i2t: failed n1");     
  if ((p[0] = (int *) malloc ((size_t) n1 * n2 * sizeof (int))) == NULL)
      failed ( (char *) "i2t: failed n2");
  for (i = 0; i < n1 - 1; i++)
    p[i + 1] = p[i] + n2;
  for (i = 0, a = p[0]; i < n1 * n2; i++)
    *a++ = 0;
  return p;
}

/*********************************/


void free_c1t(char *p){

  free(p);
}

/*********************************/


void free_i2t(int **p){

  free(p[0]);
  free(p);
}


/*********************************/
