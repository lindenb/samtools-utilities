/**
 * Author:
 *	Pierre Lindenbaum PhD
 * Contact:
 *	plindenbaum@yahoo.fr
 * Reference:
 *	http://plindenbaum.blogspot.com
 * WWW:
 *	http://plindenbaum.blogspot.com
 *	http://samtools.sourceforge.net/
 * Motivation:
 *	dump DNA genomic sequence as a CGI
 * Compilation:
 *	export SAMDIR=\$\{HOME\}/samtools-0.1.17
 *	export BUILD="\\\"toy"\\\"
 *	export GENOME_PATH="\\\"/home/lindenb/samtools-0.1.17/examples/toy.fa\\\""
 *	make bin/faidx.cgi
 */
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <errno.h>
#ifndef _NO_RAZF
#include "razf.h"
#else
#define RAZF FILE
#define razf_read(fp, buf, size) fread(buf, 1, size, fp)
#define razf_open(fn, mode) fopen(fn, mode)
#define razf_close(fp) fclose(fp)
#define razf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define razf_tell(fp) ftello(fp)
#endif

#ifndef GENOME_PATH
    #error macro GENOME_PATH was not defined (full path to a genome indexed with faidx)
#else
    #define GENOME_FAIDX GENOME_PATH ".fai"
#endif

#ifndef BUILD
   #define BUILD "undefined"
#endif

#ifndef MIN
#define MIN(a,b) (a<b?a:b)
#endif


#define NOTEMPTY(var,msg)  if(var==NULL || var[0]==0) die(msg " undefined",400)
#define CONTENT_DISPOSITION(ext)    printf("Content-Disposition: inline; filename=" BUILD "_%s_%ld_%ld.%s;\n",chrom,chromStart,chromEnd,ext)
#define HEADER(mime,ext) CONTENT_DISPOSITION(ext);printf("Content-Type:" mime ";charset=UTF-8%c%c",10,10);header_printed=1;

static const int BUFFER_LENGTH=100000;

typedef struct
	{
	int32_t line_len, line_blen;
	int64_t len;
	uint64_t offset;
	}faidx1_t;

static char* decode (char* str)
    {
    char* p=str;
    while(*p!=0)
        {
        if(*p=='+')
            {
            *p=' ';
            }
        else if(*p=='%' && *(p+1)!=0 && *(p+2)!=0)
            {
            char buffer[3]={*(p+1),*(p+2),0};
            *p=(char)strtoul(buffer,NULL,16);
            *(p+1)=0;
            p+=2;
            }
        ++p;
        }
    return str;
    }

static void die(const char*message,int status)
	{
	printf("Status: %d\n",status);
	printf("Content-Type:text/plain;charset=UTF-8%c%c",10,10);
	printf("%s.\n",message);
	fflush(stdout);
	exit(EXIT_FAILURE);
	}

static void echo(
	const char* chrom,
	long chromStart,
	long chromEnd,
	faidx1_t* index,
	int every
	)
	{
	RAZF *rz;
	int64_t pos;
	long toPrint;
	long printed=0L;	
	
	if(chromStart>= index->len) return ;
	chromEnd=MIN(index->len,chromEnd);
	if(chromStart==chromEnd) return;
	toPrint=chromEnd-chromStart;
	
	rz = razf_open(GENOME_PATH, "r");
	if (rz == NULL)
		{
		die("Cannot open " GENOME_PATH,500);
		}
	pos= 	index->offset +
		chromStart / (index->line_blen * index->line_len) +
		chromStart % index->line_blen
		;
	razf_seek(rz,pos , SEEK_SET);
	while(printed<toPrint)
		{
		long i=0;
		char buff[BUFSIZ];
		long nRead=razf_read(rz,buff,BUFSIZ);
		while(i<nRead && printed<toPrint)
			{
			if(isgraph(buff[i]))
				{
				if(every!=-1 && printed%every==0)
					{
					fputc('\n',stdout);
					}
				fputc(buff[i],stdout);
				printed++;
				}
			++i;
			}
		}
	razf_close(rz);
	}


int main(int argc,char** argv)
    {
    char* endptr=NULL;
 #ifndef TEST 
    char* method=getenv("REQUEST_METHOD");
    char* query_string=getenv("QUERY_STRING");
 #else
    char* query_string=NULL;
 #endif
    char* format=NULL;
    char* chrom=NULL;
    char* chrom_start_str=NULL;
    char* chrom_end_str=NULL;


#ifndef TEST 
    if(method==NULL || strcmp(method,"GET")!=0)
        {
        die("GET method expected",406);
        }
#endif

#ifdef TEST
   if(argc<2) return EXIT_FAILURE;
   query_string=argv[1];
#endif


    if(query_string==NULL || query_string[0]==0)
        {
        die("QUERY_STRING missing",406);
        }
    /* decode CGI string */
    const char* end=query_string+strlen(query_string);
    char* prev=query_string;
    while(prev!=end)
        {
        char* amp=strchr(prev,'&');
        char* eq=NULL;
        if(amp==NULL) amp=(char*)end;
        *amp=0;
        eq=strchr(prev,'=');
        if(eq!=NULL && eq!=prev)
            {
            char* key=prev;
            char* value;
            *eq=0;
            value=++eq;
            key=decode(key);
            value=decode(value);
            if(strcmp(key,"tid")==0 || strcmp(key,"chrom")==0)
                {
                chrom=value;
                }
            else if(strcmp(key,"start")==0)
                {
                chrom_start_str=value;
                }
            else if(strcmp(key,"end")==0)
                {
                chrom_end_str=value;
                }
            else if(strcmp(key,"fmt")==0 || strcmp(key,"format")==0)
                {
                format=value;
                }
            }
        if(amp==end) break;
        prev=++amp;
        }
    if(query_string==NULL || query_string[0]==0)
        {
        die("chrom missing",400);
        }
    NOTEMPTY(chrom,"chrom");
    NOTEMPTY(chrom_start_str,"start");
    NOTEMPTY(chrom_end_str,"end");
    errno=0;
    long chromStart=strtol(chrom_start_str,&endptr,10);
    if(chromStart<0 || errno!=0 ||*endptr!=0) die("bad value for chromStart",400);
    long chromEnd=strtol(chrom_end_str,&endptr,10);
    if(chromEnd<chromStart || errno!=0 ||*endptr!=0) die("bad value for chromEnd",400);
    
    faidx1_t index;
    index.len=-1;
    /* open index */
    FILE* faidx=fopen(GENOME_FAIDX,"r");
    if(faidx==NULL)
    	{
    	die("cannot load index " GENOME_FAIDX ".",500);
    	}
    else
    	{
    	int c;
    	int len=0;
    	size_t line_buff=BUFSIZ;
    	char* line=malloc(line_buff*sizeof(char));
    	int nameLength=strlen(chrom);
    	if(line==NULL)
    		{
    		die("Out of memory.",500);
    		}
    	while((c=fgetc(faidx))!=EOF)
    		{
    		if(c=='\n')
    			{
    			line[len]=0;
    			if(len>0 &&
    			   strncmp(chrom,line,nameLength)==0 &&
    			   line[nameLength]=='\t'
    			   )
				{
				if(sscanf(&line[nameLength+1],"%Ld\t%Ld\t%d\t%d",
					&index.len, &index.offset, &index.line_blen,&index.line_len
					)!=4)
					{
					fprintf(stderr,"#%s\n",line);
					die("Cannot read index",500);
					}
				break;
				}
    			len=0;
    			continue;	
    			}
    		if(len+2>=line_buff)
    			{
    			line_buff+=BUFSIZ;
    			line=realloc(line,line_buff*sizeof(char));
    			if(line==NULL)
		    		{
		    		die("Out of memory.",500);
		    		}
		    	
    			}
    		line[len++]=c;
    		}
    	fclose(faidx);
    	}
    if(index.len<1)
    	{
    	die("Unknown chromosome",404);
    	}
    	
    int header_printed=0;
    if(format!=NULL && strcasecmp(format,"json")==0)
        {
        HEADER("application/json","json");
        printf("{\"build\":\"%s\",\"chrom\":\"%s\",",BUILD,chrom);
        printf("\"start\":%ld,",chromStart);
        printf("\"end\":%ld,",chromEnd);
        printf("\"sequence\":\"");
        echo(chrom,chromStart,chromEnd,&index,-1);
        printf("\"}\n");
        }
    else if(format!=NULL && strcasecmp(format,"xml")==0)
        {
        HEADER("text/xml","xml");
        fputs("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n<sequence",stdout);
        printf(" build=\"%s\"",BUILD);
        printf(" chrom=\"%s\"",chrom);
        printf(" start=\"%ld\"",chromStart);
        printf(" end=\"%ld\">",chromEnd);
        echo(chrom,chromStart,chromEnd,&index,-1);
        printf("</sequence>\n");
        }
    else if(format!=NULL && strcasecmp(format,"text")==0)
        {
        HEADER("text/plain","txt");
        echo(chrom,chromStart,chromEnd,&index,-1);
        }
    else
        {
        HEADER("text/plain","fa");
        printf(">%s|%s:%ld-%ld",BUILD,chrom,chromStart,chromEnd);
        echo(chrom,chromStart,chromEnd,&index,50);
        fputc('\n',stdout);
        }
    fflush(stdout);
    return 0;
    }
