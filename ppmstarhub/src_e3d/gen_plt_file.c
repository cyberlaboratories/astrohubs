#include <stdio.h>
#include <sys/types.h>
#include <sys/file.h>
#include <math.h>

float xmin = 1.0, xmax = 0.0;

main(argc, argv)
int argc;
char *argv[];
{
  int i, iycol, nfiles, nrows, i10shift[100];
  char cOutput[256], cOption[256], coord1, cscale[256], cFiles[100][256];
  float flog10scale[100], fget_log_scale(char *cFile, int iycol);
  void gen_i10shift(int n, float *f10scale, int *i10shift);
  int get_data_file_names(char cFiles[100][256]);
  void gen_plt_file(int nfiles, char cFiles[100][256], char coord1,
                    char *cscale, int iycol, int i10shift[100]);

  // Get inputs
  if(argc < 2) { printf("usage: gen_plt_file <output>  [<option1> <option2> ...]\n"); exit(0); }
  sscanf(argv[1], "%s", cOutput);
  coord1 = cOutput[0];
  iycol = 2;
  sprintf(cscale,"linear");
   for(i=2; i<argc; i++) {
     sscanf(argv[i], "%s", cOption);
     if(strcmp(cOption,"part=sum") == 0) iycol = 7;
     if(strcmp(cOption,"logscale" ) == 0) sprintf(cscale,"%s",cOption);
     if(strcmp(cOption,"autoscale") == 0) sprintf(cscale,"%s",cOption);
  }

  // Get files, x-range & generate factor of 10 scaleings if requested & needed
  nfiles = get_data_file_names(cFiles);
  for(i=0; i<nfiles; i++) { flog10scale[i] =  fget_log_scale(cFiles[i], iycol); }
  for(i=0; i<nfiles; i++) i10shift[i] = 0;
  if(strcmp(cscale,"autoscale")==0) gen_i10shift(nfiles, flog10scale, i10shift);

  gen_plt_file(nfiles,cFiles,coord1,cscale,iycol,i10shift);
}

/****************************************************************************/
void gen_plt_file(int nfiles, char cFiles[100][256], char coord1,
                  char *cscale, int iycol, int i10shift[100])
{
  FILE *fp, *fopen();
  int i;
  char cfld[100][256], cycol[256], cyfld[256], c10[256], cXlabel[256],cuse1[256];
  void get_var_name(int i0, char cfin, char* cfile, char* cfld);
  void gen_c10(int n, char *c10);
  float xdel;
  void get_xlabel(char *cXlabel);
  void get_use1(char *cXlabel);
  void get_var_names(int nvars, char cVars[100][256]);

  sprintf(cXlabel, "set xlabel '%c'\n", coord1);
  get_xlabel(cXlabel);
  get_use1(cuse1);

  for(i=0; i<nfiles; i++) get_var_name(11,'.',cFiles[i],cfld[i]);
  get_var_names(nfiles, cfld);

  xdel = 0.02 * (xmax - xmin);
  fp = fopen("bar.plt", "w");
  // fprintf(fp,"set xlabel '%c'\n", coord1);
  fprintf(fp,"%s\n", cXlabel);
  fprintf(fp,"set ylabel '");
  for(i=0; i<nfiles; i++) fprintf(fp, " %s", cfld[i]);
  fprintf(fp,"'\n");
  if(cuse1[0] == '1') fprintf(fp,"set xrange [%f:%f]\n", xmin-xdel, xmax+xdel);
  if(strcmp(cscale,"logscale")==0) fprintf(fp, "set logscale y\n");
  fprintf(fp,"plot \\\n");
  for(i=0; i<nfiles; i++) {
    sprintf(cycol, "%d", iycol);
    sprintf(cyfld, "%s", cfld[i]);
    if(i10shift[i] > 0) {
      gen_c10(i10shift[i], c10);
      sprintf(cycol, "(%s.0*$%d)", c10, iycol);
      sprintf(cyfld, "1e%d * %s", i10shift[i], cfld[i]);
      if(i10shift[i] == 1) sprintf(cyfld, "10 * %s", cfld[i]);
      if(i10shift[i] == 2) sprintf(cyfld, "100 * %s", cfld[i]);
    }
    fprintf(fp, " '%s' using %s:%s title '%s' with linespoints", cFiles[i],cuse1,cycol,cyfld);
    if(i < nfiles-1) fprintf(fp, ", \\\n");  else fprintf(fp, "\n");
  }
  fclose(fp);
}
/****************************************************************************/
void get_xlabel(char *cXlabel)
{
  FILE    *fp;
  char    *line = NULL;
  size_t  len = 0;
  ssize_t read;

  fp = fopen("xlabel_string", "r");
  if (fp == NULL) return;

  read = getline(&line, &len, fp);
  if(read != -1) {
/*
     printf("========================================\n");
     printf("Retrieved line of length %zu :\n", read);
     printf("%s", line);
     printf("========================================\n");
*/
     sprintf(cXlabel, "%s", line);
  }
  if (line) free(line);
  fclose(fp);
  return;
}

/****************************************************************************/
void get_use1(char *cuse1)
{
  FILE    *fp;
  char    *line = NULL;
  size_t  len = 0;
  ssize_t read;
  int i, i0, i1;

  sprintf(cuse1, "1");
  fp = fopen("plot_string", "r");
  if (fp == NULL) return;

  read = getline(&line, &len, fp);
  if(read != -1) {
/*
     printf("========================================\n");
     printf("Retrieved line of length %zu :\n", read);
     printf("%s", line);
     printf("========================================\n");
*/
  i0 = -1;
  i1 = -1;
  for(i=0; i<(int)read; i++) {
    if(line[i] == '(') i0 = i;
    if(line[i] == ')') i1 = i;
  }
  if(i1 > i0  &&  i0 > 0) {
    sprintf(cuse1, "%s", &line[i0]);
    cuse1[1+i1-i0] = '\0';
  }
  }
  if (line) free(line);
  return;
}

/****************************************************************************/
void gen_c10(int n, char *c10)
{
  int i;
  c10[0] = '1';
 for(i=0; i<n; i++) c10[i+1] = '0';
 c10[n+1] = '\0';
 return;
}
/****************************************************************************/
void get_var_name(int i0, char cfin, char* cfile, char* cfld)
{
  int i=0;
  while(cfile[i0+i] != cfin  &&  i<100) { cfld[i] = cfile[i+i0]; i++; }
  cfld[i] = '\0';
  return;
}

/****************************************************************************/
void get_var_names(int nvars, char cVars[100][256])
{
  FILE    *fp;
  char    *line = NULL;
  size_t  len = 0;
  ssize_t read;
  int iVar, n;
  char ca[256], cb[256], cc[256];

  fp = fopen("ylabel_strings", "r");
  if (fp == NULL) return;

  for(iVar=0; iVar<nvars; iVar++) {
    read = getline(&line, &len, fp);
    if(read != -1) {
/*
     printf("========================================\n");
     printf("Retrieved line of length %zu :\n", read);
     printf("%s", line);
     printf("========================================\n");
*/
       sscanf(line, "%s %s %s", ca, cb, cc);
       n = strlen(cc);
       strncpy(cVars[iVar], &cc[1], n-2);
       cVars[iVar][n-1] = '\0';
    }
  }
  if (line) free(line);
  fclose(fp);
}

/****************************************************************************/
int get_data_file_names(char cFiles[100][256])
{
  FILE *fp, *fopen();
  int nfiles = 0;
  fp = fopen("plot_data_files", "r");
  while(fscanf(fp, "%s", cFiles[nfiles]) > 0) nfiles++;
  fclose(fp);
  return nfiles;
}
/****************************************************************************/
void gen_i10shift(int n, float *f10scale, int *i10shift)
{
  int i;
  float f10max;
  f10max = f10scale[0];
  for(i=1; i<n; i++) if(f10max < f10scale[i]) f10max = f10scale[i];
  for(i=0; i<n; i++) {
    i10shift[i] = 0;
    if(f10scale[i] > -8000.0) i10shift[i] = (int)(f10max+0.5 - f10scale[i]);
  }
  return;
}

/****************************************************************************/
float fget_log_scale(char *cFile, int iycol)
{
  float f10scale;
  float min, max;
  void get_minmax(char *cFile, int iycol, float *pmin, float *pmax);
  get_minmax(cFile, iycol, &min, &max);
  if(max < 0.0) max = -max;
  if(min < 0.0) min = -min;
  if(max < min) max = min;   //  = max(abs(min), abs(max))
  if(max == 0.0) f10scale = -9000.0;
  else           f10scale = log(max) / log(10.0);
  return f10scale;
}

/****************************************************************************/
#define NMAX 100000
void get_minmax(char *cFile, int iycol, float *min, float *max)
{
  int i, nrows;
  float x[NMAX], y[NMAX];
  int read_cols(char *cFileIn, int icol1, int icol2, float *x, float *y);
  nrows = read_cols(cFile, 1, iycol, x, y);
 *min = y[0];
 *max = y[0];
  if(xmax < xmin) { xmin = x[0]; xmax = x[0]; }
  for(i=0; i<nrows; i++) {
    if(*min > y[i]) *min = y[i];
    if(*max < y[i]) *max = y[i];
    if(xmin > x[i]) xmin = x[i];
    if(xmax < x[i]) xmax = x[i];
  }
}

/****************************************************************************/
int read_cols(char *cFileIn, int icol1, int icol2, float *x, float *y)
{
  FILE *fp, *fopen();
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  float d[20];
  int icolmax, i,nrows, ngot;

  if(icol1 < icol2) icolmax = icol2; else icolmax = icol1;
  fp = fopen(cFileIn, "r");
  if(fp == NULL) { printf("Could not opne %s\n", cFileIn); exit(0); }
  nrows = 0;
  while ((read = getline(&line, &len, fp)) != -1  &&  nrows < NMAX) {
    if(line[0] != '#') {
      ngot = sscanf(line, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                    &d[ 1],&d[ 2],&d[ 3],&d[ 4],&d[ 5],&d[ 6],&d[ 7],&d[ 8],&d[ 9], &d[10],
                    &d[11],&d[12],&d[13],&d[14],&d[15],&d[16],&d[17],&d[18],&d[19]);
      if(ngot < icolmax) { printf("Max column in a row of file: %d\n", ngot+1); exit(0); }
      x[nrows] = d[icol1];
      y[nrows] = d[icol2];
      nrows++;
    }
  }
  fclose(fp);
  return nrows;
}
