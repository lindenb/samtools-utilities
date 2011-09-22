/*
Motivation:
	Join a file containing some genomic positions and the
	content of a bgz file indexed with tabix
Author:
	Pierre Lindenbaum PhD
WWW:
	http://plindenbaum.blogspot.com
Contact:
	plindenbaum@yahoo.fr
Reference:
	http://plindenbaum.blogspot.com/2011/09/joining-genomic-annotations-files-with.html
Compilation:
	gcc -o jointabix -Wall -O2 -I${TABIXDIR} -L${TABIXDIR} jointabix.c  -ltabix -lz
API:
	http://samtools.sourceforge.net/tabix.shtml

*/
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include <errno.h>
#include "bgzf.h"
#include "tabix.h"

typedef struct {
	char* line;
	size_t buffer;
	size_t len;
	char** tokens;
	size_t tokens_buffer;
	size_t n_tokens;
	char delim; 
	gzFile in;
	char ignore;
	int chromCol;
	int startCol;
	int endCol;
	int shift;
	tabix_t *t;
	} JoinTabix;


static char* readline(JoinTabix* app)
	{
	int c;
	app->len=0;
	if(gzeof(app->in)) return NULL;
	
	
	while((c=gzgetc(app->in))!=EOF && c!='\n')
		{
		if(app->len+2 >= app->buffer)
			{
			app->buffer+=BUFSIZ;
			if((app->line=(char*)realloc(app->line,sizeof(char)*app->buffer))==NULL)
				{
				fprintf(stderr,"Cannot realloc %d bytes\n",app->buffer);
				exit(EXIT_FAILURE);
				}
			}

		app->line[app->len++]=c;	
		}
	app->line[app->len]=0;
	return app->line;
	}

static void pushToken(JoinTabix* app,char* ptr)
	{
	if(app->n_tokens+1 > app->tokens_buffer)
		{
		app->tokens_buffer+=50;
		if((app->tokens=(char**)realloc(app->tokens,sizeof(char*)*app->tokens_buffer))==NULL)
			{
			fprintf(stderr,"Cannot realloc %d bytes\n",sizeof(char*)*app->tokens_buffer);
			exit(EXIT_FAILURE);
			}
		}
	
	app->tokens[app->n_tokens++]=ptr;	
	}

static void splitLine(JoinTabix* app)
	{
	size_t i,j;

	app->n_tokens=0;
	pushToken(app,app->line);
	for(i=0, j=1;i< app->len;++i)
		{
		if(app->line[i]!=app->delim) continue;
		app->line[i]=0;
		pushToken(app,&app->line[i+1]);
		}
	}

static void printTokens(FILE* out,JoinTabix* app)
	{
	size_t i;

	for(i=0;i< app->n_tokens;++i)
		{
		if(i>0) fputc(app->delim,out);
		fputs(app->tokens[i],out);
		}
	}
static int parseIntGE0(const char* s)
	{
	char* p2;
	int n;

	errno=0;
	n=(int)strtol(s,&p2,10);
	if(*p2!=0 || errno!=0 || n<0)
		{
		return -1;
		}
	return n;
	}
	
	
static int parseInt(const char* s)
	{
	int n=parseIntGE0(s);
	if(n<0)
		{
		fprintf(stderr,"Bad Integer \"%s\".\n",s);
		return EXIT_FAILURE;
		}
	return n;
	}

static int parseInt1(const char* s)
	{
	int n=parseInt(s);
	if(n<1)
		{
		fprintf(stderr,"Expected a number greater than 0 . Got \"%s\".\n",s);
		return EXIT_FAILURE;
		}
	return n-1;/* convert to 0-based */
	}


static int join(JoinTabix* app)
 	{
 	while(readline(app)!=NULL)
 		{
 		int tid=0;
 		int found=0;
 		int chromStart=0;
 		int chromEnd=0;
 		if(app->line[0]==0) continue;
 		
 		if(app->line[0]==app->ignore)
 			{
 			fputs(app->line,stdout);
 			fputc('\n',stdout);
 			continue;
 			}
 		splitLine(app);

 		
 		if(app->chromCol>=app->n_tokens ||
 		   app->startCol>=app->n_tokens ||
 		   app->endCol>=app->n_tokens ||
 		   (tid=ti_get_tid(app->t->idx, app->tokens[app->chromCol])) <0 ||
 		   (chromStart=parseIntGE0(app->tokens[app->startCol])) <0 ||
 		   (chromEnd=parseIntGE0(app->tokens[app->endCol])) <0
 		   )
 			{
 			fprintf(stderr,"Found column missing or bad position or unknown chromosome in : ");
 			printTokens(stderr,app);
 			fputc('\n',stderr);
 			}
 		else
 			{
 			const char *s;
 			int len;
 			if(chromStart==chromEnd) ++chromEnd;
 			chromStart+= app->shift;
 			chromEnd+= app->shift;
 			ti_iter_t iter = ti_queryi(app->t, tid, chromStart, chromEnd);
			while ((s = ti_read(app->t, iter, &len)) != 0)
				{
				printTokens(stdout,app);
				fputc(app->delim, stdout);
				fputs(s, stdout);
				fputc('\n', stdout);
				++found;
				}
			ti_iter_destroy(iter);
			}
		if(found==0)
			{
			fputs("##boum\t",stdout);
 			printTokens(stdout,app);
 			fputc('\n',stdout);
			}
 		}
 	return 0;
 	}


int main(int argc, char *argv[])
  {
  char* tabixfile=NULL;
  JoinTabix param;
  int optind=1;
  memset((void*)&param,0,sizeof(JoinTabix));
  param.delim='\t';
  param.chromCol=0;
  param.startCol=1;
  param.endCol=1;
  param.ignore='#';
  /* loop over the arguments */
  while(optind<argc)
	    {
	    if(strcmp(argv[optind],"-h")==0)
		    {
		    fprintf(stdout, "Author: Pierre Lindenbaum PHD. 2011.\n");
		    fprintf(stdout, "Last compilation:%s %s\n",__DATE__,__TIME__);
		    fprintf(stdout, "Usage: %s (options) {stdin|file|gzfiles}:\n",argv[0]);
		    fprintf(stdout, "  -d <char> column delimiter. default: TAB\n");
		    fprintf(stdout, "  -c <int> chromosome column (%d).\n",param.chromCol+1);
		    fprintf(stdout, "  -s <int> start column (%d).\n",param.startCol+1);
		    fprintf(stdout, "  -e <int> end column (%d).\n",param.endCol+1);
		    fprintf(stdout, "  -i <char> ignore lines starting with (\'%c\').\n",param.ignore);
		    fprintf(stdout, "  -t <filename> tabix file (required).\n");
		    fprintf(stdout, "  +1 add 1 to the genomic coodinates.\n");
		    fprintf(stdout, "  -1 remove 1 to the genomic coodinates.\n");
		    return EXIT_SUCCESS;
		    }
	    else if(strcmp(argv[optind],"-f")==0 && optind+1< argc)
	    	{
	    	tabixfile=argv[++optind];
	    	}
	    else if(strcmp(argv[optind],"-1")==0)
	    	{
	    	param.shift=-1;
	    	}
	    else if(strcmp(argv[optind],"+1")==0)
	    	{
	    	param.shift=1;
	    	}
	    else if(strcmp(argv[optind],"-d")==0 && optind+1< argc)
		    {
		    if(strlen(argv[optind+1])!=1)
		    	{
		    	fprintf(stderr,"Expected only one char for the delimiter.\n");
		    	return EXIT_FAILURE;
		    	}
		    param.delim=argv[++optind][0];
		    }
	    else if(strcmp(argv[optind],"-c")==0 && optind+1< argc)
	    	{
	    	param.chromCol=parseInt1(argv[++optind]);
	    	}
	    else if(strcmp(argv[optind],"-s")==0 && optind+1< argc)
	    	{
	    	param.startCol=parseInt1(argv[++optind]);
	    	}
	    else if(strcmp(argv[optind],"-e")==0 && optind+1< argc)
	    	{
	    	param.endCol=parseInt1(argv[++optind]);
	    	}
	    else if(strcmp(argv[optind],"-i")==0 && optind+1< argc)
	    	{
	    	if(strlen(argv[optind+1])!=1)
		    	{
		    	fprintf(stderr,"Expected only one char for the delimiter.\n");
		    	return EXIT_FAILURE;
		    	}
	    	param.ignore=argv[++optind][0];
	    	}
	    else if(strcmp(argv[optind],"--")==0)
		    {
		    optind++;
		    break;
		    }
	    else if(argv[optind][0]=='-')
		    {
		    fprintf(stderr,"Unnown option: %s\n",argv[optind]);
		    return EXIT_FAILURE;
		    }
	    else
		    {
		    break;
		    }
	    ++optind;
	    }
  if(tabixfile==NULL)
	{
	fprintf(stderr,"Error: undefined tabix file.\n");
	return EXIT_FAILURE;
	}
	    
  if(param.chromCol==param.startCol)
	{
	fprintf(stderr,"Error: col(chrom)==col(start).\n");
	return EXIT_FAILURE;
	}

  if(param.chromCol==param.endCol)
	{
	fprintf(stderr,"Error: col(chrom)==col(end).\n");
	return EXIT_FAILURE;
	}

  if ((param.t = ti_open(tabixfile, 0)) == 0)
	{
	fprintf(stderr, "Cannot open tabix file \"%s\" %s.\n",tabixfile,strerror(errno));
	return EXIT_FAILURE;
	}
  if (ti_lazy_index_load(param.t) < 0)
  	 {
         fprintf(stderr, "Cannot open index for file \"%s\".\n",tabixfile);
	 return EXIT_FAILURE;
         }

  if(optind==argc)
      {
      param.in=gzdopen(fileno(stdin),"r");
      if(param.in==NULL)
      	{
        fprintf(stderr,"Cannot read from stdin\n");
        return EXIT_FAILURE;
      	}
      join(&param);
      }
  /* loop over the files */
  while(optind<argc)
	{
	char* filename=argv[optind++];
	errno=0;
	param.in=gzopen(filename,"r");
	if(param.in==NULL)
      		{
        	fprintf(stderr,"Cannot read %s (%s).\n",filename,strerror(errno));
        	return EXIT_FAILURE;
      		}
      	join(&param);
      	gzclose(param.in); 
	}
  /* we're done */
  ti_close(param.t);
  free(param.line);
  free(param.tokens);
  return EXIT_SUCCESS;
  }
