/**
 * Author:
 *	Pierre Lindenbaum PhD
 * Contact:
 *	plindenbaum@yahoo.fr
 * WWW:
 *	http://plindenbaum.blogspot.com
 *	http://samtools.sourceforge.net/
 *	http://samtools.sourceforge.net/sam-exam.shtml
 * Motivation:
 *	creates a WIG file from a BAM file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "sam.h"

typedef struct parameter_t{
	FILE* out;
	int prev_tid;
	int prev_pos;
	int beg;
	int end;
	samfile_t *in;
} Param,*ParamPtr;

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
	bam_plbuf_t *buf = (bam_plbuf_t*)data;
	bam_plbuf_push(b, buf);
	return 0;
}
// callback for bam_plbuf_init()
static int  scan_all_genome_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
	{
	ParamPtr param = (ParamPtr)data;
	if ((int)pos >= param->beg && (int)pos < param->end)
		{
		if(param->prev_tid!=tid || param->prev_pos+1!=(int)pos)
			{
			fprintf(param->out,"name=%s pos=%d\n", param->in->header->target_name[tid],pos);
			}
	
		fprintf(param->out,"%d\n",n);
		param->prev_tid=(int)tid;
		param->prev_pos=(int)pos;
		}
	
	return 0;
	}

static void usage()
	{
	fprintf(stdout, "Author: Pierre Lindenbaum PHD. 2011.\n");
	fprintf(stdout, "Last compilation:%s %s\n",__DATE__,__TIME__);
	fprintf(stdout, "Usage: bam2wig (options) <aln.bam>\n");
	fprintf(stdout, "Options:\n");
	}
	
int main(int argc, char *argv[])
	{
	int optind=1;
	Param parameter;
	while(optind < argc)
		{
		if(strcmp(argv[optind],"-h")==0)
		        {
		        usage();
		        return EXIT_FAILURE;
		        }
		else if(strcmp(argv[optind],"-g")==0 && optind+1<argc)
		        {
		      
		        }
		else if(strcmp(argv[optind],"--")==0)
		        {
		        ++optind;
		        break;
		        }
		else if(argv[optind][0]=='-')
		        {
		        fprintf(stderr,"%s: unknown option '%s'\n",argv[0],argv[optind]);
		        exit(EXIT_FAILURE);
		        }
		else
		        {
		        break;
		        }
		++optind;
		}
	
	
	
	
	parameter.out=stdout;
	parameter.prev_tid=-1;
	parameter.prev_pos=-1;
	parameter.beg = 0;
	parameter.end = INT_MAX;
	
	parameter.in = samopen(argv[optind], "rb", 0);
	if (parameter.in == 0)
		{
		fprintf(stderr, "Cannot open BAM file \"%s\".\n", argv[optind]);
		return EXIT_FAILURE;
		}
	

	
	if (argc == 2)
		{
		sampileup(parameter.in, -1, scan_all_genome_func, &parameter);
		}	
	else {
		int ref;
		bam_index_t *idx;
		bam_plbuf_t *buf;
		idx = bam_index_load(argv[1]); // load BAM index
		if (idx == 0) {
			fprintf(stderr, "BAM indexing file is not available.\n");
			return 1;
		}
		bam_parse_region(parameter.in->header, argv[2], &ref,
		                 &parameter.beg, &parameter.end); // parse the region
		if (ref < 0) {
			fprintf(stderr, "Invalid region %s\n", argv[2]);
			return 1;
		}
		buf = bam_plbuf_init( scan_all_genome_func, &parameter); // initialize pileup
		bam_fetch(parameter.in->x.bam, idx, ref, parameter.beg, parameter.end, buf, fetch_func);
		bam_plbuf_push(0, buf); // finalize pileup
		bam_index_destroy(idx);
		bam_plbuf_destroy(buf);
	}
	samclose(parameter.in);
	return 0;
	}

