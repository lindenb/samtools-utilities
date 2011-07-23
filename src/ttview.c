/**
 * Author:
 *	Pierre Lindenbaum PhD
 * Orginal source code:
 *	bam_tview.c in http://samtools.sourceforge.net/
 *	Authors: Heng Li, Bob Handsaker, Jue Ruan, Colin Hercus, Petr Danecek
 * Contact:
 *	plindenbaum@yahoo.fr
 * Reference:
 *	http://plindenbaum.blogspot.com/2011/07/text-alignment-viewer-using-samtools.html
 * WWW:
 *	http://plindenbaum.blogspot.com
 *	http://samtools.sourceforge.net/
 * Compilation:
 *	make ##(generate samtools)
 *	gcc -o bamttview -g -Wall -O2 -DSTANDALONE_VERSION  -I. -Lbcftools  bam_ttview.c  bam2bcf.o   errmod.o  bam_color.o libbam.a -lbcf  -lm -lz
 * Motivation:
 *	Text alignment viewer using the samtools API
 */
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <errno.h>
#include <limits.h>
#include "bam.h"
#include "faidx.h"
#include "bam2bcf.h"

char bam_aux_getCEi(bam1_t *b, int i);
char bam_aux_getCSi(bam1_t *b, int i);
char bam_aux_getCQi(bam1_t *b, int i);

#define TV_MIN_ALNROW 2
#define TV_MAX_GOTO  40
#define TV_LOW_MAPQ  10

#define TV_COLOR_MAPQ   0
#define TV_COLOR_BASEQ  1
#define TV_COLOR_NUCL   2
#define TV_COLOR_COL    3
#define TV_COLOR_COLQ   4

#define TV_BASE_NUCL 0
#define TV_BASE_COLOR_SPACE 1

/* I put the characters in a structure. If one day I want to add colors... */
typedef struct cpixel_t
	{
	char c;
	}CPixel;

typedef struct
	{
	/* number of columns */
	int mcol;
	/* number of lines created so far */
	int nLines;
	/* the array(y,x) */
	CPixel** screen;

	bam_index_t *idx;
	bam_lplbuf_t *lplbuf;
	bam_header_t *header;
	bamFile fp;
	int curr_tid, left_pos;
	faidx_t *fai;
	bcf_callaux_t *bca;

	int ccol, last_pos, row_shift, base_for, color_for, is_dot, l_ref, ins, no_skip, show_name;
	char *ref;
} ttview_t;

#define WHERE fprintf(stderr,"[DEBUG] %d\n",__LINE__)


static CPixel* getchxy(ttview_t* t,int y,int x)
	{
	int i=0;
	assert(y>=0);
	if(y<0 || x<0 || x>= t->mcol ) return NULL;
	while(t->nLines<=y)
		{
		t->screen=(CPixel**)realloc(t->screen,sizeof(CPixel*)*(t->nLines+1));
		if(t->screen==NULL)
			{
			fputs("Out of memory\n",stderr);
			exit(EXIT_FAILURE);
			}
		t->screen[t->nLines]=malloc(t->mcol*sizeof(CPixel));
		if(t->screen[t->nLines]==NULL)
			{
			fputs("Out of memory\n",stderr);
			exit(EXIT_FAILURE);
			}
		for(i=0;i< t->mcol;++i)
			{
			t->screen[t->nLines][i].c=' ';
			}
		t->nLines++;
		}

	
	return &(t->screen[y][x]);
	}

static void putchxy(ttview_t* t,int y,int x,char c)
	{
	CPixel* pixel=NULL;
	pixel=getchxy(t,y,x);
	if(pixel==NULL)
		{
		return;
		}
	pixel->c=c;
	}

static void printfyx(ttview_t* t,int y,int x,char* fmt,...)
	{
	va_list ap;    
	int i=0;
	char* buffer=NULL;

	buffer=malloc(t->mcol+1);

	if(buffer==NULL)
		{
		fputs("Out of memory\n",stderr);
		exit(EXIT_FAILURE);
		}
	memset(buffer,'\0',sizeof(char)*(t->mcol+1));

	va_start(ap,fmt);
	vsnprintf(buffer,t->mcol,fmt,ap);
	va_end(ap);
	while(x+i< t->mcol && buffer[i]!=0)
		{
		putchxy(t,y,x+i,buffer[i]);
		++i;
		}
	free(buffer);
	}

static void dump(ttview_t* t)
	{
	int y;
	int x;
	for(y=0;y< t->nLines;++y)
		{
		for(x=0;x< t->mcol;++x)
			{
			fputc(t->screen[y][x].c,stdout);
			}
		fputc('\n',stdout);
		}
	}


static int ttv_pl_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
	{
	extern unsigned char bam_nt16_table[256];
	ttview_t *tv = (ttview_t*)data;
	int i, j, c, rb, max_ins = 0;
	uint32_t call = 0;
	if (pos < tv->left_pos || tv->ccol > tv->mcol) return 0; // out of screen
	// print reference
	rb = (tv->ref && pos - tv->left_pos < tv->l_ref)? tv->ref[pos - tv->left_pos] : 'N';
	for (i = tv->last_pos + 1; i < pos; ++i) {
		if (i%10 == 0 && tv->mcol - tv->ccol >= 10) printfyx(tv,0, tv->ccol, "%-d", i+1);
		c = tv->ref? tv->ref[i - tv->left_pos] : 'N';
		putchxy(tv,1, tv->ccol++, c);
	}
	if (pos%10 == 0 && tv->mcol - tv->ccol >= 10) printfyx(tv,0, tv->ccol, "%-d", pos+1);
	{ // call consensus
		bcf_callret1_t bcr;
		int qsum[4], a1, a2, tmp;
		double p[3], prior = 30;
		bcf_call_glfgen(n, pl, bam_nt16_table[rb], tv->bca, &bcr);
		for (i = 0; i < 4; ++i) qsum[i] = bcr.qsum[i]<<2 | i;
		for (i = 1; i < 4; ++i) // insertion sort
			for (j = i; j > 0 && qsum[j] > qsum[j-1]; --j)
				tmp = qsum[j], qsum[j] = qsum[j-1], qsum[j-1] = tmp;
		a1 = qsum[0]&3; a2 = qsum[1]&3;
		p[0] = bcr.p[a1*5+a1]; p[1] = bcr.p[a1*5+a2] + prior; p[2] = bcr.p[a2*5+a2];
		if ("ACGT"[a1] != toupper(rb)) p[0] += prior + 3;
		if ("ACGT"[a2] != toupper(rb)) p[2] += prior + 3;
		if (p[0] < p[1] && p[0] < p[2]) call = (1<<a1)<<16 | (int)((p[1]<p[2]?p[1]:p[2]) - p[0] + .499);
		else if (p[2] < p[1] && p[2] < p[0]) call = (1<<a2)<<16 | (int)((p[0]<p[1]?p[0]:p[1]) - p[2] + .499);
		else call = (1<<a1|1<<a2)<<16 | (int)((p[0]<p[2]?p[0]:p[2]) - p[1] + .499);
	}
	
	c = ",ACMGRSVTWYHKDBN"[call>>16&0xf];
	i = (call&0xffff)/10+1;
	if (i > 4) i = 4;
	if (c == toupper(rb)) c = '.';
	putchxy(tv,2, tv->ccol, c);
	if(tv->ins) {
		// calculate maximum insert
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = pl + i;
			if (p->indel > 0 && max_ins < p->indel) max_ins = p->indel;
		}
	}
	// core loop
	for (j = 0; j <= max_ins; ++j) {
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = pl + i;
			int row = TV_MIN_ALNROW + p->level - tv->row_shift;
			if (j == 0) {
				if (!p->is_del) {
					if (tv->base_for == TV_BASE_COLOR_SPACE && 
							(c = bam_aux_getCSi(p->b, p->qpos))) {
						c = bam_aux_getCSi(p->b, p->qpos);
						// assume that if we found one color, we will be able to get the color error
						if (tv->is_dot && '-' == bam_aux_getCEi(p->b, p->qpos)) c = bam1_strand(p->b)? ',' : '.';
					} else {
						if (tv->show_name) {
							char *name = bam1_qname(p->b);
							c = (p->qpos + 1 >= p->b->core.l_qname)? ' ' : name[p->qpos];
						} else {
							c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
							if (tv->is_dot && toupper(c) == toupper(rb)) c = bam1_strand(p->b)? ',' : '.';
						}
					}
				} else c = p->is_refskip? (bam1_strand(p->b)? '<' : '>') : '*';
			} else { // padding
				if (j > p->indel) c = '*';
				else { // insertion
					if (tv->base_for ==  TV_BASE_NUCL) {
						if (tv->show_name) {
							char *name = bam1_qname(p->b);
							c = (p->qpos + j + 1 >= p->b->core.l_qname)? ' ' : name[p->qpos + j];
						} else {
							c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
							if (j == 0 && tv->is_dot && toupper(c) == toupper(rb)) c = bam1_strand(p->b)? ',' : '.';
						}
					} else {
						c = bam_aux_getCSi(p->b, p->qpos + j);
						if (tv->is_dot && '-' == bam_aux_getCEi(p->b, p->qpos + j)) c = bam1_strand(p->b)? ',' : '.';
					}
				}
			}
			if (row > TV_MIN_ALNROW
				/* && row < tv->mrow */
				)
				{
				int x;
				
				
				if (tv->color_for == TV_COLOR_BASEQ)
					{
					x = bam1_qual(p->b)[p->qpos]/10 + 1;
					if (x > 4) x = 4;	
					} 
				else if (tv->color_for == TV_COLOR_MAPQ)
					{
					x = p->b->core.qual/10 + 1;
					if (x > 4) x = 4;
					}	
				else if (tv->color_for == TV_COLOR_NUCL)
					{
					x = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)] + 5;
					}	
				else if(tv->color_for == TV_COLOR_COL)
					{
					x = 0;
					switch(bam_aux_getCSi(p->b, p->qpos))
						{
						case '0': x = 0; break;
						case '1': x = 1; break;
						case '2': x = 2; break;
						case '3': x = 3; break;
						case '4': x = 4; break;
						default: x = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)]; break;
						}
					x+=5;
					}	
				else if(tv->color_for == TV_COLOR_COLQ)
					{
					x = bam_aux_getCQi(p->b, p->qpos);
					if(0 == x) x = bam1_qual(p->b)[p->qpos];
					x = x/10 + 1;
					if (x > 4) x = 4;
					}
				
				putchxy(tv,row, tv->ccol, bam1_strand(p->b)? tolower(c) : toupper(c));
				
			}
		}
		c = j? '*' : rb;
		if (c == '*') {
			
			
			putchxy(tv,1, tv->ccol++, c);
			
		} else putchxy(tv,1, tv->ccol++, c);
	}
	tv->last_pos = pos;
	return 0;
}

static ttview_t *ttv_init(const char *fn, const char *fn_fa)
	{
	ttview_t *tv = (ttview_t*)calloc(1, sizeof(ttview_t));
	if(tv==NULL) return NULL;
	tv->screen=NULL;
	tv->nLines=0;
	tv->is_dot = 1;
	tv->fp = bam_open(fn, "r");
	bgzf_set_cache_size(tv->fp, 8 * 1024 *1024);
	assert(tv->fp);
	tv->header = bam_header_read(tv->fp);
	tv->idx = bam_index_load(fn);
	if (tv->idx == 0) exit(1);
	tv->lplbuf = bam_lplbuf_init(ttv_pl_func, tv);
	if (fn_fa) tv->fai = fai_load(fn_fa);
	tv->bca = bcf_call_init(0.83, 13);
	tv->ins = 1;

	
	tv->mcol = 80;

	
	tv->color_for = TV_COLOR_MAPQ;
	
	return tv;
	}
static void ttv_clear(ttview_t *tv)
	{
	int i=0;
	for(i=0;i< tv->nLines;++i)
		{
		free(tv->screen[i]);
		}
	
	free(tv->screen);
	tv->screen=NULL;
	tv->nLines=0;
	}
static void ttv_destroy(ttview_t *tv)
	{
	ttv_clear(tv);
	bam_lplbuf_destroy(tv->lplbuf);
	bcf_call_destroy(tv->bca);
	bam_index_destroy(tv->idx);
	if (tv->fai) fai_destroy(tv->fai);
	free(tv->ref);
	bam_header_destroy(tv->header);
	bam_close(tv->fp);
	free(tv);
	}

static int ttv_fetch_func(const bam1_t *b, void *data)
{
	ttview_t *tv = (ttview_t*)data;
	if (tv->no_skip) {
		uint32_t *cigar = bam1_cigar(b); // this is cheating...
		int i;
		for (i = 0; i <b->core.n_cigar; ++i) {
			if ((cigar[i]&0xf) == BAM_CREF_SKIP)
				cigar[i] = cigar[i]>>4<<4 | BAM_CDEL;
		}
	}
	bam_lplbuf_push(b, tv->lplbuf);
	return 0;
}

static int ttv_draw_aln(ttview_t *tv, int tid, int pos)
	{
	
	// reset
	ttv_clear(tv);
	tv->curr_tid = tid;
	tv->left_pos = pos;
	tv->last_pos = tv->left_pos - 1;
	tv->ccol = 0;
	// print ref and consensus
	if (tv->fai) {
		char *str;
		if (tv->ref) free(tv->ref);
		str = (char*)calloc(strlen(tv->header->target_name[tv->curr_tid]) + 30, 1);
		sprintf(str, "%s:%d-%d", tv->header->target_name[tv->curr_tid], tv->left_pos + 1, tv->left_pos + tv->mcol);
		tv->ref = fai_fetch(tv->fai, str, &tv->l_ref);
		free(str);
	}
	// draw aln
	bam_lplbuf_reset(tv->lplbuf);
	bam_fetch(tv->fp, tv->idx, tv->curr_tid, tv->left_pos, tv->left_pos + tv->mcol, tv, ttv_fetch_func);
	bam_lplbuf_push(0, tv->lplbuf);

	while (tv->ccol < tv->mcol)
		{
		int pos = tv->last_pos + 1;
		if (pos%10 == 0 && tv->mcol - tv->ccol >= 10) printfyx(tv,0, tv->ccol, "%-d", pos+1);
		putchxy(tv,1, tv->ccol++, (tv->ref && pos < tv->l_ref)? tv->ref[pos - tv->left_pos] : 'N');
		++tv->last_pos;
		}
	return 0;
	}


static void usage()
	{
	 fprintf(stderr,"Pierre Lindenbaum PHD. 2010.\nOrginal code: from the samtools package http://samtools.sourceforge.net . \n");
	fprintf(stderr, "Usage: bamtk ttview (options) <aln.bam> [ref.fasta]\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -g <region>\n");
	fprintf(stderr, "  -f <filename> reads a list regions ( '-' for stdin)\n");
	fprintf(stderr, "  -d toggle dot view\n");
	//fprintf(stderr, "  -s ref skip\n");
	fprintf(stderr, "  -r toggle read name\n");
	fprintf(stderr, "  -N nt view\n");
	fprintf(stderr, "  -C cs view\n");
	fprintf(stderr, "  -i insertions\n");
	fprintf(stderr, "  -X <int> number of columns\n");
	fprintf(stderr, "  -T <positive int> shift all positions by <T>.\n");
	}

int bam_ttview_main(int argc, char *argv[])
	{
	int optind=1;
	ttview_t *tv;
	char* region=NULL;
	char* filename=NULL;
	int shift=0;
	int base_for=TV_BASE_NUCL, is_dot=1, ins=1, no_skip=0, show_name=0;
	int columns=80;
	 while(optind < argc)
		{
		if(strcmp(argv[optind],"-h")==0)
		        {
		        usage();
		        exit(EXIT_FAILURE);
		        }
		else if(strcmp(argv[optind],"-g")==0 && optind+1<argc)
		        {
		        region=argv[++optind];
		        }
		else if(strcmp(argv[optind],"-X")==0 && optind+1<argc)
		        {
		        columns=atoi(argv[++optind]);
		        if(columns<=10) columns=10;
		        }
		else if(strcmp(argv[optind],"-T")==0 && optind+1<argc)
		        {
		        shift=atoi(argv[++optind]);
		        if(shift<=0) shift=0;
		        }
		else if(strcmp(argv[optind],"-f")==0 && optind+1<argc)
		        {
		        filename=argv[++optind];
		        }
		else if(strcmp(argv[optind],"-d")==0 )
		        {
		        is_dot=!is_dot;
		        }
		else if(strcmp(argv[optind],"-N")==0 )
		        {
		        base_for=TV_BASE_NUCL;
		        }
		else if(strcmp(argv[optind],"-C")==0 )
			{
			base_for=TV_BASE_COLOR_SPACE;
			} 
		else if(strcmp(argv[optind],"-i")==0 )
			{
			ins=!ins;
			}
		else if(strcmp(argv[optind],"-r")==0 )
			{
			show_name=!show_name;
			} 
		else if(strcmp(argv[optind],"-s")==0 )
			{
			no_skip=!no_skip;
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
	if(!(optind+2==argc || optind+1==argc))
		{
		usage();
		exit(EXIT_FAILURE);
		}
	
	tv = ttv_init(argv[optind], (argc != optind+2)? 0 : argv[optind+1]);
	if(tv==NULL) return -1;
	tv->base_for=base_for;
	tv->is_dot=is_dot;
	tv->ins=ins;
	tv->no_skip=no_skip;
	tv->show_name=show_name;
	tv->mcol=columns;
	if(region!=NULL)
		{
		int tid = -1, beg,end;
		bam_parse_region(tv->header, region, &tid, &beg, &end);
		
		if (tid < 0 || tid>=tv->header->n_targets)
			{
			fprintf(stderr,"Bad region %s\n",region);
			}
		else
			{
			ttv_draw_aln(tv, tid,  (beg-shift<0?0:beg-shift));
			dump(tv);
			}
		}
	else if(filename!=NULL)
		{
		int tid = -1, beg,end;
		char line[BUFSIZ];
		FILE* in=stdin;
		if(strcmp(filename,"-")!=0)
			{
			in=fopen(filename,"r");
			if(in==NULL)
				{
				fprintf(stderr,"Cannot open %s %s\n",filename,strerror(errno));
				exit(EXIT_FAILURE);
				}
			}
		while(fgets(line,BUFSIZ,in)!=NULL)
			{
			if(line[0]=='\n' || line[0]=='#') continue;
			
			bam_parse_region(tv->header, line, &tid, &beg, &end);
			
			if (tid < 0 || tid>=tv->header->n_targets)
				{
				fprintf(stderr,"Bad region %s\n",region);
				}
			else
				{
				ttv_draw_aln(tv, tid, (beg-shift<0?0:beg-shift));
				fprintf(stdout,"\n\n> %s\n",line);
				dump(tv);
				fputc('\n',stdout);
				}
			}
		if(strcmp(filename,"-")!=0)
			{
			fclose(in);
			}
		}
	else
		{
		ttv_draw_aln(tv,0,0);
		dump(tv);
		}
	ttv_destroy(tv);
	return 0;
	}

#ifdef STANDALONE_VERSION
int main(int argc,char** argv)
	{
	return bam_ttview_main(argc,argv);
	}
#endif
