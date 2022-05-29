/*
 * FIXME: Not until now (Jan 29, 2022) I realised that i ignored
 * the collisions that happens because the particles are deviated
 * by previous collisions and happens outside the window.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define Box(i,j,k) boxes[(i*n_boxes+j)*n_boxes+k]

#define get_part_idx(pt) (pt-particles)

typedef double vec[3];

typedef struct Particle* part;

struct Particle{
    vec x;
    vec v;
    part n,p;
    double last;
    int flag;
    int i[3];
};

typedef struct Event{
    enum {part_t, wall_t, box_t, sys_t} type;
    double time, reg_time;
    union{
        struct{
            part p[2];
            vec  v[2];
        }p;
        struct{
            part p;
            int  w;
        }w;
        struct{
            part p;
            int  b;
        }b;
        void(*s)(void);
    }e;
} event;

event*  events;
event*  last_event;
part    particles;
part*   boxes;
int     n_particles;
int     n_boxes;
int     n_boxes_tot;
int     n_events;
int     max_events;
int     flag;
int     debug;
double  current;
double  included;
double  radius;
double  fix_t;
double  length;
double  max_dt;
double  min_dt;
double  box_size;
double  max_velocity;
double  final_t;

void new_part(){
    particles = malloc(sizeof(struct Particle)*n_particles);
    for(int i=0; i<n_particles; i++){
        for(int j=0; j<3; j++){
            particles[i].x[j] = radius+(double)rand()/RAND_MAX*(length-2*radius);
            particles[i].v[j] = (double)rand()/RAND_MAX*max_velocity*2-max_velocity;
        }
        particles[i].p = particles[i].n = particles + i;
        particles[i].last = -1;
        particles[i].flag = 0;
    }
}

void add_part(part new, int i, int j, int k){ 
    part box  = Box(i,j,k);
    if(box==NULL)
        Box(i,j,k) = new;
    else{
        new->n    = box->n;
        new->n->p = new;
        new->p    = box;
        box->n    = new;
        new->i[0] = i;
        new->i[1] = j;
        new->i[2] = k;
    }
}

void rm_part(part x){
    if(Box(x->i[0],x->i[1],x->i[2])==x)
        if(x->n == x)
            Box(x->i[0],x->i[1],x->i[2]) = NULL;
        else
            Box(x->i[0],x->i[1],x->i[2]) = x->n;
    x->p->n = x->n;
    x->n->p = x->p;
    x->n = x->p = x;
}

void part_into_box(){
    printf("entering part_into_box... ");
    for(int n=0; n<n_particles; n++){
        int i=particles[n].x[0] / box_size;
        int j=particles[n].x[1] / box_size;
        int k=particles[n].x[2] / box_size;
        add_part(particles+n, i,j,k);
    }
    printf("done\n");
}

void add_event(){
    if (n_events == max_events){
        printf("heap overflow: too many events\n");
        exit(1);
    }
    event new = events[n_events];
    int i = n_events++;
    while(i!=0 && events[(i-1)/2].time > new.time){
        events[i] = events[(i-1)/2];
        i = (i-1)/2;
    }
    events[i] = new;
}

void pop_event(){
    printf("pop,%d...", n_events);
    event poped = events[0];
    int i=0;
    while(2*i+2 < n_events){
        int min = events[2*i+1].time < events[2*i+2].time? 2*i+1: 2*i+2;
        events[i] = events[min];
        i = min;
    }
    if(2*i+2==n_events)
        events[i] = events[n_events-1];
    else{
        event new = events[n_events-1];
        while(i!=0 && events[(i-1)/2].time > new.time){
            events[i] = events[(i-1)/2];
            i = (i-1)/2;
        }
        events[i] = new;
    }
    events[--n_events] = poped;
    printf("done\n");
}

void try_part(part a, part b){
    vec relx = {(a->x[0]-b->x[0])/2,
                (a->x[1]-b->x[1])/2,
                (a->x[2]-b->x[2])/2};
    vec relv = {(a->v[0]-b->v[0])/2,
                (a->v[1]-b->v[1])/2,
                (a->v[2]-b->v[2])/2};
    double cost;
    if((cost=(relx[0]*relv[0]+relx[1]*relv[1]+relx[2]*relv[2]))>0)
        return;
    double relx2 = relx[0]*relx[0]+relx[1]*relx[1]+relx[2]*relx[2];
    double relv2 = relv[0]*relv[0]+relv[1]*relv[1]+relv[2]*relv[2];
    double relx_n = sqrt(relx2);
    double relv_n = sqrt(relv2);
    cost /= (relx_n * relv_n);
    double delta = relx2*(cost*cost-1)+radius*radius;
    if (delta<0) return;
    double d = -relx_n*cost-sqrt(delta);
    double dt = d / relv_n;
    if (dt > max_dt || dt + current < included) return;
    for(int i=0; i<3; i++) relx[i] += dt * relv[i];
    double s = relx[0]*relx[0]+relx[1]*relx[1]+relx[2]*relx[2];
    s = 2*(relx[0]*relv[0]+relx[1]*relv[1]+relx[2]*relv[2])/s;
    for(int i=0; i<3; i++) relv[i] = s*relx[i] - relv[i];
    event* e = events+n_events;
    e->type        = part_t;
    e->reg_time    = current;
    e->time        = current + dt;
    e->e.p.p[0]    = a;
    e->e.p.p[1]    = b;
    e->e.p.v[0][0] = (a->v[0]+b->v[0])/2 + relv[0];
    e->e.p.v[0][1] = (a->v[1]+b->v[1])/2 + relv[1];
    e->e.p.v[0][2] = (a->v[2]+b->v[2])/2 + relv[2];
    e->e.p.v[1][0] = (a->v[0]+b->v[0])/2 - relv[0];
    e->e.p.v[1][1] = (a->v[1]+b->v[1])/2 - relv[1];
    e->e.p.v[1][2] = (a->v[2]+b->v[2])/2 - relv[2];
    add_event();
}

void try_wall(part a){
    double dt=1000;
    int w;
    for(int i=0; i<3; i++){
        double _dt = (radius+(length-2*radius)*(a->v[i] > 0) - a->x[i]) / a->v[i];
        if(_dt < dt){
            dt = _dt;
            w  = i;
        }
    }
    if(dt > max_dt || current + dt < included) return;
    event*e = events+n_events;
    e->type  = wall_t;
    e->reg_time = current;
    e->time  = current + dt;
    e->e.w.p = a;
    e->e.w.w = w;
    add_event();
}

void try_box (part a){
    double dt = 1000;
    int b;
    for(int i=0; i<3; i++){
        int foo = a->v[i] > 0;
        double _dt = ( (a->i[i] + foo) * box_size - a->x[i]) / a->v[i];
        if(_dt<dt){
            dt = _dt;
            b  = foo?1+i:-1-i;
        }
    }
    if(dt>max_dt || dt+current<=included
        || b==-1 && a->x[0]<box_size
        || b==1  && a->x[0]>length-box_size
        || b==-2 && a->x[1]<box_size
        || b==2  && a->x[1]>length-box_size
        || b==-3 && a->x[2]<box_size
        || b==3  && a->x[2]>length-box_size) return;
    event*e = events+n_events;
    e->type  = box_t;
    e->reg_time = current;
    e->time  = dt + current;
    e->e.b.p = a;
    e->e.b.b = b;
    add_event();
}

void try(){
    printf("trying...");
    if(last_event==NULL || last_event->type==sys_t){
	for(int i=0; i<n_boxes_tot; i++) if(boxes[i]!=NULL){
	    part a = boxes[i];
	    part _a = a;
	    do{
		a->flag = flag;
		part b = a->n;
		while(b!=_a){
		    try_part(a,b);
		    b = b->n;
		}
		try_wall(a);
		try_box(a);
		a=a->n;
	    }while(a != _a);
	}
    }else if(last_event->type==wall_t || last_event->type==box_t){
	part a = last_event->type==wall_t
	    ? last_event->e.w.p
	    : last_event->e.b.p;
	try_wall(a);
	try_box(a);
	part b = Box(a->i[0],a->i[1],a->i[2]);
	part _b = b;
	do{
	    if(b!=a) try_part(a,b);
	    b = b->n;
	}while(b!=_b);
    }else{ // part_t
	part a = last_event->e.p.p[0];
	part b = last_event->e.p.p[1];
	try_wall(a); try_box(a);
	try_wall(b); try_box(b);
	try_part(a,b);
	part c = Box(a->i[0],a->i[1],a->i[2]);
	part _c = c;
	do{
	    if(c!=a) try_part(a,c);
	    if(c!=b) try_part(b,c);
	    c = c->n;
	}while(c!=_c);
    }
    included = current + max_dt;
    flag = 1-flag;
    printf("done\n");
}
 
const double tol = 1e-6;

int check_mutual(part a){
    part box = Box(a->i[0],a->i[1],a->i[2]);
    part cur = box;
    do {
        if ( cur == a ) return 1;
        cur = cur -> n;
    } while(cur != box);
    return 0;
}

int check_box(){
    int errcnt=0;
    int tot_checked=0;
    for(int i=0; i<n_boxes; i++)
        for(int j=0; j<n_boxes; j++)
            for(int k=0; k<n_boxes; k++)
                if(Box(i,j,k)!=NULL){
                    part start = Box(i,j,k);
                    part cur = start;
                    do{
                        if(cur->x[0] > length+tol || cur->x[0] < 0-tol ||
                           cur->x[1] > length+tol || cur->x[1] < 0-tol ||
                           cur->x[2] > length+tol || cur->x[2] < 0-tol )
                                {printf("out of bound,%ld(%lf,%lf,%lf),mutual:%d\n",cur-particles,cur->x[0],cur->x[1],cur->x[2],check_mutual(cur));
                                 exit(1);}
                        if(cur->x[0] > (i+1)*box_size+tol ||
                           cur->x[0] <     i*box_size-tol ||
                           cur->x[1] > (j+1)*box_size+tol ||
                           cur->x[1] <     j*box_size-tol ||
                           cur->x[2] > (k+1)*box_size+tol ||
                           cur->x[2] <     k*box_size-tol){
                            errcnt ++;
#if 0
                            printf("box err: part %ld x=(%lf,%lf,%lf), box=(%d,%d,%d)\n",
                                    cur - particles, cur->x[0], cur->x[1], cur->x[2], i, j, k);
#endif
                        }
                        cur = cur->n;
			tot_checked ++;
                    }while(cur!=start);
                }
    printf("checked=%d,tot=%d\n", tot_checked, n_particles);
    return errcnt;
}

void fix_box(){
    printf("fix box...");
    for(int n=0; n<n_particles; n++){
        if( particles[n].x[0] > (1+particles[n].i[0])*box_size+tol || 
            particles[n].x[0] <    particles[n].i[0] *box_size-tol ||
            particles[n].x[1] > (1+particles[n].i[1])*box_size+tol || 
            particles[n].x[1] <    particles[n].i[1] *box_size-tol ||
            particles[n].x[2] > (1+particles[n].i[2])*box_size+tol || 
            particles[n].x[2] <    particles[n].i[2] *box_size-tol ) {
                printf("fix: %d\n", n);
                rm_part(particles+n);
                add_part(particles+n,(int)(particles[n].x[0]/box_size),
                                     (int)(particles[n].x[1]/box_size),
                                     (int)(particles[n].x[2]/box_size));
        }
    }
    printf("donefix:%d\n",check_box());
}

int main(int argc, char *argv[]){
    n_particles =1000;
    max_events  = 10000;
    n_boxes     = 5;
    n_boxes_tot = n_boxes*n_boxes*n_boxes;
    radius      = .5;
    length      = 10;
    max_velocity = 1;
    final_t     = 10;
    max_dt      = .1;
    min_dt      = 1e-10;
    debug       = 0;
    fix_t       = 1;

    int opt;
    while((opt = getopt(argc, argv, "n:e:b:r:l:v:T:t:f:dhs")) != -1)
        switch (opt) {
            case 'n': n_particles = atoi(optarg); break;
            case 'e': max_events  = atoi(optarg); break;
            case 'b': n_boxes     = atoi(optarg); break;
            case 'r': radius      = atof(optarg); break;
            case 'l': length      = atof(optarg); break;
            case 'v': max_velocity = atof(optarg); break;
            case 'T': final_t     = atof(optarg); break;
            case 't': max_dt      = atof(optarg); break;
            case 'f': fix_t       = atof(optarg); break;
            case 'd': debug       = 1;            break;
	    case 's': srand(time(NULL));	  break;
            case 'h':
	    default : printf("\nrigid: boltzmann simulation. Usage:\n"
			     "  -n number of particles\n"
                             "  -e maximum number of events\n"
                             "  -b number of boxes\n"
                             "  -r radius of each particle\n"
                             "  -l length of the space\n"
                             "  -v maximum velocity when initializing\n"
                             "  -T final time\n"
                             "  -t maximum step\n"
                             "  -f interval of fixing the boxes\n"
                             "  -d debug\n"
                             "  -s use randome seed\n"
                             "  -h show this message\n\n"); exit(1);
        }

    box_size    = length/n_boxes;
    n_events    = 0;
    flag        = 1;
    current     = 0;
    included    = 0;
    last_event	= NULL;
    
    new_part();

    events      = malloc(sizeof(event)*max_events);
    boxes       = malloc(sizeof(part)*n_boxes_tot);
    double dt   = 0;
    int step    = 0, _step;
    int part_cnt= 0, wall_cnt=0, box_cnt=0;
    for(int i=0; i<n_boxes_tot; i++) boxes[i] = NULL;

    part_into_box();

    for(double _t=fix_t; _t<final_t; _t+=fix_t){
        events[n_events].time = _t;
        events[n_events].type = sys_t;
        events[n_events].e.s  = &fix_box;
        add_event();
    }

    event* e=NULL;
    while(current<final_t){

	int _check_box = check_box();
	printf("t = %lf, check_box() = %d\n", current, _check_box);
	//if(_check_box>0)exit(1);

        _step = current / max_dt;
        if(_step>step){
            step=_step;
            printf("[ %2.10f ] part: %d, wall: %d, box: %d\n", 
                    current,part_cnt,wall_cnt,box_cnt);
        }

        try();

        printf("n_event:%d\n", n_events);
        do{
            if( n_events==0 ){
                e = NULL;
                break;
            }
            pop_event();
            e = events + n_events;
        }while(e->type==part_t && 
                  ((e->reg_time < e->e.p.p[0]->last || e->time < e->e.p.p[1]->last))
               || e->type==wall_t && e->reg_time < e->e.w.p->last
               || e->type==box_t  && e->reg_time < e->e.b.p->last
               || isnan(e->time)
	       || (e->time-current) < min_dt
               || e->time < 0);

	last_event = e;

        if (e!=NULL) 
            dt = e->time - current;
        else
            dt = max_dt;
	printf("dt=%16.14lf\n", dt);
        for(int i=0; i<n_particles; i++)
            for(int j=0; j<3; j++)
                particles[i].x[j] += particles[i].v[j] * dt;

        current += dt;
 
        if (debug) {
            printf("%.12f ", current);
            if(e==NULL)
                printf("no event\n");
            else if (e->type==part_t)
                printf("part %ld %ld x1=(%lf,%lf,%lf) v1=(%lf,%lf,%lf) x2=(%lf,%lf,%lf) v2=(%lf,%lf,%lf)\n",
                        get_part_idx(e->e.p.p[0]), get_part_idx(e->e.p.p[1]),
                        e->e.p.p[0]->x[0],e->e.p.p[0]->x[1],e->e.p.p[0]->x[2],
                        e->e.p.p[0]->v[0],e->e.p.p[0]->v[1],e->e.p.p[0]->v[2],
                        e->e.p.p[1]->x[0],e->e.p.p[1]->x[1],e->e.p.p[1]->x[2],
                        e->e.p.p[1]->v[0],e->e.p.p[1]->v[1],e->e.p.p[1]->v[2]);
            else if (e->type==box_t)
                printf("box %ld x=(%lf,%lf,%lf) v=(%lf,%lf,%lf) b=%d",
                        get_part_idx(e->e.b.p),
                        e->e.b.p->x[0],e->e.b.p->x[1],e->e.b.p->x[2],
                        e->e.b.p->v[0],e->e.b.p->v[1],e->e.b.p->v[2], e->e.b.b);
            else if (e->type==wall_t)
                printf("wall %ld x=(%lf,%lf,%lf) v=(%lf,%lf,%lf) w=%d\n",
                        get_part_idx(e->e.w.p),
                        e->e.w.p->x[0],e->e.w.p->x[1],e->e.w.p->x[2],
                        e->e.w.p->v[0],e->e.w.p->v[1],e->e.w.p->v[2], e->e.w.w);
            else
                printf("sys \n");
            printf("errcnt=%d\n",check_box());
        }

        if(e->type==part_t){
            for(int j=0; j<3; j++){
                e->e.p.p[0]->v[j] = e->e.p.v[0][j];
                e->e.p.p[1]->v[j] = e->e.p.v[1][j];
            }
            e->e.p.p[0]->last = e->e.p.p[1]->last = current;
            part_cnt++;
        }else if(e->type==wall_t){
            e->e.w.p->v[ e->e.w.w ] *= -1;
            e->e.w.p->last = current;
            wall_cnt++;
        }else if(e->type==box_t){
            int i = e->e.b.p->i[0];
            int j = e->e.b.p->i[1];
            int k = e->e.b.p->i[2];
            printf("(%d,%d,%d)->",i,j,k);
            i += (e->e.b.b==1) ? 1 : ( (e->e.b.b==-1)? -1 : 0 );
            j += (e->e.b.b==2) ? 1 : ( (e->e.b.b==-2)? -1 : 0 );
            k += (e->e.b.b==3) ? 1 : ( (e->e.b.b==-3)? -1 : 0 );
            printf("(%d,%d,%d)\n",i,j,k);
            if( !( (i<0)||(i==n_boxes)||(j<0)||(j==n_boxes)||(k<0)||(k==n_boxes) ) ){
                rm_part(e->e.b.p);
                add_part(e->e.b.p, i,j,k);
                box_cnt++;
            }
        }else if(e->type==sys_t)
            e->e.s();
    }
    return 0;
}
