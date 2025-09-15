
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<ctype.h>
#include <unistd.h>
/**
 *@file: CalVSP.c  
 *                                                                                                
 *CalVSP: a program for analyzing the molecular surface areas, volumes, and polar surface areas            
 * @author YuZhu Li
 * @email	drugyuzhu@163.com
 * @date 2025-09-15
 * 
 * @copyright Copyright (c) 2025   Shanghai Jiao Tong University && China State Institute of Pharmaceutical Industry
 * This project is released under the MIT License.
 * See LICENSE file for details.                                                                                                                
*/








struct min_max {float xmin;float xmax;float ymin;float ymax;float zmin;float zmax;};
struct grid {float x;float y; float z; int state;short int psa;};
struct data_index{int x;int y;int z;int wsp;};
struct atoms {float xa; float ya; float za;int ele;float rvdw;};
struct atom_info{ int falg; int id; float vdwr; float mw; char atoms[5];};


struct grid *grid_data; 
struct atoms *atom_data;
struct min_max  get_min_max_mol2(char *name);
struct min_max  get_min_max(struct atoms *atom_data,int atom_number);
struct result{ float surf; float psa; float vol;};
int cal_dis(struct atoms *atom_data,int atom_number,float water_radius,float x,float y,float z);
int cal_adjacent(int number, struct grid *grid_data,int grid_number,int x_number,int y_number,int z_number);
int cal_adjacent1(int number, struct grid *grid_data,int grid_number,int x_number,int y_number,int z_number);
float  cal_grid_dis(struct grid grid_dataa,struct grid grid_datab);
int cal_psa(struct grid grid_data,struct atoms *atom_data,int atom_number);
int file_type(char *file_name);
int get_atom_info(struct atom_info *atom_info_data);
int get_data(struct atoms **atom_data,struct atom_info *atom_info_data,int count, float x,float y,float z,char *ele_s);
int get_data_from_file(char *file_name,struct atom_info *atom_info_data,struct atoms **atom_data,int *p_atom_number,int *xyz_file_number);
void trace_falef(struct atoms *atom_data_trace,struct atoms *atom_data, int atom_number,int file_number);
int calculate(struct atoms *atom_data, struct result *resultdata,int if_trace_file,int atom_number,float water_flag,float water_radius,int file_number,float *surface_data,float *psa_data,float *volume_data,int out_data);
int get_xyz_atom_number(char *name);
int get_xyz_file_number(char *name);
int is_obabel_installed();
int cal_confab(char *file_name);
int cal_conformer(char *file_name);
int cal_obabel_mmf(char *file_name);
int get_mopac_dataf(char *file_name,char *out_name);
int cal_mopac_mmf(char *file_name,int multiplicity,int charge, float permittivity,char *mopac_path );
int get_energy(char *file_name,float bloz);
int get_dist_data(char *file_name,float *distribution, int input_filenumber );
int ask_if_bloz();
float get_dis();
void add_prefix(char *newnamef, char *name, char *add);
void add_suffix(char *newnamef, char *name, char *add);
void set_name(int cal_type,char *input_name,char *seted_name);
void chang_suffix2xyz(char *name);
int is_number(const char *str); 




int main(int argc, char *argv[])
{
	
	

	int cal_type=0;
    char *input_file = NULL,*cwater_flag = NULL,*cwater_radius = NULL,*cmopac_path_flag =NULL;
	int c_genf=0,c_method_mm=0 ; 
    
    char mopacphath[100];  // Set default MOPAC output file
    int opt,if_trace_file=0,genc_method=0,mm_method=0,mopac_multi=1,mopac_charge=0,out_data=0;  
	float water_flag=0,dielectric=0.0; // dielectric[Dielectric constant of solution]
//Short option description:
//- i <file>: Required input file 
//- d [file]: Optional output file (default output.data) 
//- s<number>: Optional SASA calculation input 
// -g <number>: Generate conformation through obabel 0 Monte Carlo method 1 confab
// -m <number> : Energy minimization method 0 MMF94 1 MOPAC 3 MMF94+MOPAC

  
    while ((opt = getopt(argc, argv, "i:s::dgm")) != -1) {  
        switch (opt) {
            case 'i':
                input_file = optarg;
                break;
            case 'd':
                out_data=1;  
                break;
        case 's': {
            char *arg = optarg;
            if (arg == NULL) {
                if (optind < argc && is_number(argv[optind])) {
                    arg = argv[optind];
                    optind++;  
                } else {
                    water_flag = 1.4;
                    break;
                }
            }

            if (*arg == '=') arg++;

            char *endptr;
            double s_value = strtod(arg, &endptr);
            if (endptr == arg || *endptr != '\0') {
                fprintf(stderr, "Error: -s The parameter must be a floating-point number\n");
                return 1;
            }
            water_flag = s_value;
            break;
        } 
       		case 'g':
            	c_genf = 1; 
            break;
        	case 'm':
            	c_method_mm= 1; 
            break;
            default:
                fprintf(stderr, "Usage: %s -i input file [-d output file] [-s SASA value] [-r The radius of water molecules]\n", argv[0]);
                return 1;
        }
    }


	printf("water flag ridus %f\n",water_flag);
    // Check required parameters
    if (input_file == NULL) {
        fprintf(stderr, "Error: Input file must be specified (-i)\n");
        return 1;
    }

/**************************Optional parameters********************************/


// generation conforatiom
	char buffer_command[300];
	if(c_genf!=0){
		
		printf("There are two methods for generation conformation 1 Monte Carlo method 2 confab\n");
		printf("Recommend using [2]confab, which runs slowly but generates fewer wrong stractures.\n");
		printf("Please enter the selected method:\n");		
		fgets(buffer_command, sizeof(buffer_command), stdin);
		sscanf(buffer_command, "%d", &genc_method);	
		if(genc_method==1) printf("You chooed method is£º[1] Monte Carlo method \n");	
		else if(genc_method==2) printf("You chooed method is£º[2] confab method \n");
		else{
			printf("You chooed method is£º[%], undefined exit\n");
			return -1;
		}		
	}
// Energy minimization method 0 MMF94 1 MOPAC 3 MMF94+MOPAC
	if(c_method_mm!=0){
		printf("There are three methods for generation conformation 1:MMF94 2:MOPAC 3:MMF94+MOPAC\n");
		printf("Please enter the selected method:\n");		
		fgets(buffer_command, sizeof(buffer_command), stdin);
		sscanf(buffer_command, "%d", &mm_method);
		if((mm_method==2)||(mm_method==3)){
			if(mm_method==2)
			printf("You chooed method for generation conformation is£º[%d]:MOPAC \n",mm_method);
			if(mm_method==3)
			printf("You chooed method for generation conformation is£º[%d]:MMF94+MOPAC \n",mm_method);
			printf("Please enter the multiplicity for MOPAC calculation This method only supports 1-5:\n");
			fgets(buffer_command, sizeof(buffer_command), stdin);
			sscanf(buffer_command, "%d", &mopac_multi);						
			printf("Please enter the charge for MOPAC calculation:\n");	
			fgets(buffer_command, sizeof(buffer_command), stdin);
			sscanf(buffer_command, "%d", &mopac_charge);
			printf("You entered charge is£º%d \n",mopac_charge);
			printf("Please enter the Dielectric constant of solution for MOPAC calculation\n");	
			printf("vacuum ¦Å=0.0 water ¦Å=78.3553, N-Octanol ¦Å=9.8629, DMSO£º¦Å=46.826, Heptane: ¦Å=1.9113 ...\n");	
			fgets(buffer_command, sizeof(buffer_command), stdin);
			sscanf(buffer_command, "%f", &dielectric);
			printf("You entered Dielectric constant is£º%f \n",dielectric);	
			
			printf("Please enter the path of MOPAC or Press the Enter key\n");		
			fgets(mopacphath, sizeof(buffer_command), stdin);
			int i=0;
			while(mopacphath[i]!='\n') i++;
			mopacphath[i]='\0';
			int ifpath=0;
			i=0;
			while(mopacphath[i]!='\0'){
				if((mopacphath[i]=='\\')||(mopacphath[i]==':')||(mopacphath[i]=='/'))
				ifpath++;
				i++;
			}
			if(ifpath==0){
				if(system(mopacphath)==1) ifpath=1;
			}
			if(ifpath>=1)
				printf("You entered path of mopac is£º%s \n",mopacphath);
			else{
				printf("You pressd the Enter key or  entered path of mopac is£º%s is wrong use [\\MOPAC\\MOPAC2016.exe]\n",mopacphath);	
				#ifdef _WIN32
    			// Windows  
    				strcat(mopacphath,".\\MOPAC\\mopac2016.exe");
				#else
    			// Unix/Linux/macOS  
    				strcat(mopacphath,"./MOPAC/mopac2016.exe");
				#endif
			//	exit(-1);
			}										
		}else if(mm_method==1){
			printf("You chooed method for generation conformation is£º[%d]:MMF94 \n",mm_method);
		}
		else{
			printf("You chooed method for generation conformation is£º[%d], undefined exit\n",mm_method);
			return -1;
		}
	}
/**
**calculation 
**/
	char newname[100], newname1[100], newname2[100], newname3[100], newname4[100], newname5[100];
	int if_blot=0;
	float sel_bloz=0.0;
	
	if((c_genf==0)&&(c_method_mm==0)){
		cal_type=1;	
	}else if((c_genf!=0)&&(c_method_mm==0)){
		cal_type=2;
		if(genc_method==1){
	 		cal_conformer(input_file);			
		}
		if(genc_method==2){	
			cal_confab(input_file);	
		}		
	}else if((c_genf==0)&&(c_method_mm!=0)){
		if(mm_method==1){	
			cal_obabel_mmf(input_file);
			if_blot=ask_if_bloz();
			if(if_blot==0) cal_type=4;
			else if(if_blot>0){
				cal_type=3;
				sel_bloz= get_dis();
				add_prefix(newname, input_file, "mm_");
				get_energy(newname,sel_bloz);	
			}		 
		}		
		if(mm_method==2){
			if_blot=ask_if_bloz();
			if(if_blot>0){
				sel_bloz= get_dis();
				get_energy(input_file,sel_bloz);
				add_suffix(newname,input_file, "_choosed.xyz");
				cal_mopac_mmf(newname,mopac_multi,mopac_charge, dielectric,mopacphath );
				if_blot=ask_if_bloz();
				if(if_blot>0){
					sel_bloz= get_dis();
					add_prefix(newname1, newname, "mopac_");
					get_energy(newname1,sel_bloz);
					cal_type=19;					
				}else if(if_blot==0){
					cal_type=20;
				}
									
			}else if(if_blot==0){
				cal_mopac_mmf(input_file,mopac_multi,mopac_charge, dielectric,mopacphath );
				if_blot=ask_if_bloz();
				if(if_blot==0) cal_type=6;
				else if(if_blot>0){
					cal_type=5;
					sel_bloz= get_dis();
					add_prefix(newname, input_file, "mopac_");
					get_energy(newname,sel_bloz);	
				}				
			}
		}			

		if(mm_method==3){
			cal_obabel_mmf(input_file);	
			add_prefix(newname, input_file, "mm_");
			if_blot=ask_if_bloz();
			if(if_blot==0){
				cal_mopac_mmf(newname,mopac_multi,mopac_charge, dielectric,mopacphath ); 	
				if_blot=ask_if_bloz();
				if(if_blot==0){
					cal_type=10;
				}
				else if(if_blot>0){
					cal_type=9;
					sel_bloz= get_dis();
					add_prefix(newname1, newname, "mopac_");
					get_energy(newname1,sel_bloz);
					
				}
			}
			else if(if_blot>0){
				sel_bloz= get_dis();
				get_energy(newname,sel_bloz);
				add_suffix(newname1, newname, "_choosed.xyz");
				cal_mopac_mmf(newname1,mopac_multi,mopac_charge, dielectric,mopacphath ); 
				if_blot=ask_if_bloz();
				if(if_blot==0){
					cal_type=8;
				}else if(if_blot>0){
					cal_type=7;
					sel_bloz= get_dis();
					add_prefix(newname2, newname1, "mopac_");
					get_energy(newname2,sel_bloz);
					
				}
			}
		}	
	}else if((c_genf!=0)&&(c_method_mm!=0)){
		if(genc_method==1){
	 		cal_conformer(input_file);			
		}
		if(genc_method==2){	
			cal_confab(input_file);	
		}
		add_prefix(newname, input_file, "cos_");
		chang_suffix2xyz(newname);
		if(mm_method==1){	
			cal_obabel_mmf(newname);
			if_blot=ask_if_bloz();
			if(if_blot==0) cal_type=12;
			else if(if_blot>0){
				cal_type=11;
				sel_bloz= get_dis();
				add_prefix(newname1, newname, "mm_");
				get_energy(newname1,sel_bloz);	
			}		 
		}		
		
		if(mm_method==2){	
			cal_mopac_mmf(newname,mopac_multi,mopac_charge, dielectric,mopacphath );
			if_blot=ask_if_bloz();
			if(if_blot==0) cal_type=14;
			else if(if_blot>0){
				cal_type=13;
				sel_bloz= get_dis();
				add_prefix(newname1, newname, "mopac_");
				get_energy(newname1,sel_bloz);	
			}		 
		}				
		
		if(mm_method==3){
			cal_obabel_mmf(newname);	
			add_prefix(newname1, newname, "mm_");
			if_blot=ask_if_bloz();
			if(if_blot==0){
				cal_mopac_mmf(newname1,mopac_multi,mopac_charge, dielectric,mopacphath );
				if_blot=ask_if_bloz();
				if(if_blot==0){
					cal_type=18;
				}
				else if(if_blot>0){
					cal_type=17;
					sel_bloz= get_dis();
					add_prefix(newname2, newname1, "mopac_");
					get_energy(newname2,sel_bloz);	
				}
			}
			else if(if_blot>0){
				sel_bloz= get_dis();
				get_energy(newname1,sel_bloz);
				add_suffix(newname2, newname1, "_choosed.xyz");
				cal_mopac_mmf(newname2,mopac_multi,mopac_charge, dielectric,mopacphath ); 
				if_blot=ask_if_bloz();
				if(if_blot==0){
					cal_type=16;
				}else if(if_blot>0){
					cal_type=15;
					sel_bloz= get_dis();
					add_prefix(newname3, newname2, "mopac_");
					get_energy(newname3,sel_bloz);
					
				}
			}
		}			
		
		
				
	}
//	printf("\n  calculation file type is %d \n\n",cal_type);	
/****************************get calculated file name***************/
char setname[100];
set_name(cal_type,input_file,setname);
printf("\n  calculation file is %s \n\n",setname);
/**
**Read Numbe Element vdW Radius Atomic Weight and other data of atoms from the ele_data.txt file
**/		
int atom_number,i, file_number=0;
struct atom_info *atom_info_data;
atom_info_data=(struct atom_info *)malloc(200*sizeof(struct atom_info));
get_atom_info(atom_info_data);

/**
**Read info from the input file
**/	

get_data_from_file(setname,atom_info_data, &atom_data,&atom_number,&file_number);
if(file_number>1) if_trace_file=1;
else if_trace_file=0;

struct result *frame_result;
frame_result=(struct result  *)malloc(file_number*sizeof(struct result  ));	

float surface,volume,psa;

if((if_trace_file)&&(file_type(setname)==4)){
	struct atoms *atom_data_split;
	atom_data_split=(struct atoms *)malloc(atom_number*sizeof(struct atoms ));
	for(i=1;i<=file_number;i++){
		trace_falef(atom_data,atom_data_split , atom_number, i);
		calculate(atom_data_split,frame_result,if_trace_file,atom_number,water_flag,water_flag,i,&surface,&psa,&volume,out_data);
	
	}
	free(atom_data_split);

}
else calculate(atom_data,frame_result,if_trace_file,atom_number,water_flag,water_flag,0,&surface,&psa,&volume,out_data);

if((if_trace_file)>0){
	float asuf=0,apsa=0,avol=0;
	for(i=0;i<=file_number-1;i++){
		asuf=asuf+frame_result[i].surf;
		apsa=apsa+frame_result[i].psa;
		avol=avol+frame_result[i].vol;
	}
	asuf=asuf/file_number;
	apsa=apsa/file_number;
	avol=avol/file_number;
	printf("aeverage surface: %f\naeverage psa: %f\naeverage vol: %f\n",asuf,apsa,avol);
}
free(atom_data);
free(frame_result);






return 0;
}
















float  cal_grid_dis(struct grid grid_dataa,struct grid grid_datab)
{
	return sqrt(  (grid_dataa.x-grid_datab.x)*(grid_dataa.x-grid_datab.x) + (grid_dataa.y-grid_datab.y)*(grid_dataa.y-grid_datab.y)  + (grid_dataa.z-grid_datab.z)*(grid_dataa.z-grid_datab.z));
}


/*  van der Waals radii cited from Dalton Trans., 2013, 42, 8617¨C8636*/

int cal_dis(struct atoms *atom_data,int atom_number,float water_radius,float x,float y,float z)
{

int i,fit=0 ,index_min;
float distance,min_dis;
min_dis=sqrt(  (atom_data[0].xa-x)*(atom_data[0].xa-x) + (atom_data[0].ya-y)*(atom_data[0].ya-y)  +  (atom_data[0].za- z)*(atom_data[0].za-z));
index_min=0;
for(i=0;i<=atom_number-1;i++)
{	distance=sqrt(  (atom_data[i].xa-x)*(atom_data[i].xa-x) + (atom_data[i].ya-y)*(atom_data[i].ya-y)  +  (atom_data[i].za- z)*(atom_data[i].za-z));
	if(min_dis>distance){
		min_dis=distance;
		index_min=i;	
	} 
	if(distance<=(atom_data[i].rvdw+water_radius)) fit++;
}

if(fit==0){
	if(min_dis<=(atom_data[index_min].rvdw+water_radius+2.5))
	return 0;
	else return 100;
} 
if(fit==1) return 1;
if(fit>1) return  2;
else return -1;
}


struct min_max  get_min_max(struct atoms *atom_data,int atom_number)
{
	struct min_max minmax;
	int i;
	float xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
	xmin=atom_data[0].xa;xmax=atom_data[0].xa;
	ymin=atom_data[0].ya;ymax=atom_data[0].ya;
	zmin=atom_data[0].za;zmax=atom_data[0].za;
	for(i=0;i<=atom_number-1;i++)
	{
		xmin= (atom_data[i].xa < xmin) ? atom_data[i].xa : xmin;
		xmax= (atom_data[i].xa > xmax) ? atom_data[i].xa : xmax;
		ymin= (atom_data[i].ya < ymin) ? atom_data[i].ya : ymin;
		ymax= (atom_data[i].ya > ymax) ? atom_data[i].ya : ymax;
		zmin= (atom_data[i].za < zmin) ? atom_data[i].za : zmin;
		zmax= (atom_data[i].za > zmax) ? atom_data[i].za : zmax;		
	}
	minmax.xmin=xmin;
	minmax.xmax=xmax;
	minmax.ymin=ymin;
	minmax.ymax=ymax;
	minmax.zmin=zmin;
	minmax.zmax=zmax;
	return minmax;
}




int cal_adjacent(int number, struct grid *grid_data,int grid_number,int x_number,int y_number,int z_number)
{
	int xadd1,xsub1,yadd1,ysub1,zadd1,zsub1,adja=0,count=0,gen_number1=0,gen_number2=0;
	if(grid_data[number].state==0)
	{
		
		//  x-x
		xadd1=number+(y_number*z_number);
		xsub1=number-(y_number*z_number);
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}
		


		//y-y
		yadd1=number+z_number;
		ysub1=number-z_number;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0)&&(grid_data[yadd1].state!=100)&&(grid_data[ysub1].state!=100))
			count++;
		}
		
	

		
		//z-z
		zadd1=number+1;
		zsub1=number-1;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0)&&(grid_data[zadd1].state!=100)&&(grid_data[zsub1].state!=100))
			count++;
		}	
		
		//x-1,y+1,z---------x+1,y-1,z
		xadd1=number+(y_number*z_number)-z_number;
		xsub1=number-(y_number*z_number)+z_number;		
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}
		
		
				
		// x-1,y-1,z---------x+1,y+1,z
		xadd1=number+(y_number*z_number)+z_number;
		xsub1=number-(y_number*z_number)-z_number;		
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}	
		
		
			
		// x+1,y+1,z+1--------	x-1,y-1,z-1	
		xadd1=number+(y_number*z_number)+z_number+1;
		xsub1=number-(y_number*z_number)-z_number+1;
		
	
					
		if(    ((xadd1>0)&&(xadd1<=grid_number-1)) &&
		       ((xsub1>0)&&(xsub1<=grid_number-1))
		  )
		{
		
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				
				count++;
			}
			
		}		
	

		// x,y+1,z+1--------	x,y-1,z-1	
		xadd1=number+z_number+1;
		xsub1=number-z_number-1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}		

		// x-1,y+1,z+1--------	x+1,y-1,z-1	
		xadd1=number+(y_number*z_number)-z_number-1;
		xsub1=number-(y_number*z_number)+z_number+1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}		

		// x-1,y,z+1--------	x+1,y,z-1	
		xadd1=number+(y_number*z_number)-1;
		xsub1=number-(y_number*z_number)+1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}	

		// x-1,y-1,z+1--------	x+1,y+1,z-1	
		xadd1=number+(y_number*z_number)+z_number-1;
		xsub1=number-(y_number*z_number)-z_number+1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}

		// x,y-1,z+1--------	x,y+1,z-1	
		xadd1=number+z_number-1;
		xsub1=number-z_number+1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}


		// x+1,y-1,z+1--------	x-1,y+1,z-1	
		xadd1=number+(y_number*z_number)-z_number+1;
		xsub1=number-(y_number*z_number)+z_number-1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}

		// x+1,y,z+1--------	x-1,y,z-1	
		xadd1=number+(y_number*z_number)+1;
		xsub1=number-(y_number*z_number)-1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}


		if(count>=1) return 3;
		else return 1;	
	}
	else return 0;

}


int cal_adjacent1(int number, struct grid *grid_data,int grid_number,int x_number,int y_number,int z_number)
{
	int xadd1,xsub1,yadd1,ysub1,zadd1,zsub1,adja=0,count=0,gen_number1=0,gen_number2=0;
	if(grid_data[number].state==0)
	{
		//  x-x
		xadd1=number+(y_number*z_number);
		xsub1=number-(y_number*z_number);
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}
		
		xadd1=number+(y_number*z_number)*1;
		xsub1=number-(y_number*z_number)*2;
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}
			
		
		xadd1=number+(y_number*z_number)*2;
		xsub1=number-(y_number*z_number)*1;
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}

		//y-y
		yadd1=number+z_number;
		ysub1=number-z_number;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0) )
			count++;
		}
		
		yadd1=number+z_number*2;
		ysub1=number-z_number*2;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0) )
			count++;
		}		

		yadd1=number+z_number*1;
		ysub1=number-z_number*2;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0) )
			count++;
		}		
				
		yadd1=number+z_number*2;
		ysub1=number-z_number*1;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0) )
			count++;
		}				
		//z-z
		zadd1=number+1;
		zsub1=number-1;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0) )
			count++;
		}	
		zadd1=number+2;
		zsub1=number-2;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0) )
			count++;
		}	
		zadd1=number+1;
		zsub1=number-2;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0) )
			count++;
		}	
		
						
		zadd1=number+2;
		zsub1=number-1;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0) )
			count++;
		}
		if(count>=1)
			return 3;
		else return 0;	
	}
	else return 0;
}





int cal_psa(struct grid grid_data,struct atoms *atom_data,int atom_number)
{
	int i,minatom,j,ortho_atom;
	float distance,mindis;

	distance=sqrt(  (atom_data[0].xa-grid_data.x)*(atom_data[0].xa-grid_data.x) + (atom_data[0].ya-grid_data.y)*(atom_data[0].ya-grid_data.y)  +  (atom_data[0].za- grid_data.z)*(atom_data[0].za-grid_data.z));	
	minatom=0;
	mindis=distance;
	for(i=0;i<= atom_number-1;i++)
	{
		distance=sqrt(  (atom_data[i].xa-grid_data.x)*(atom_data[i].xa-grid_data.x) + (atom_data[i].ya-grid_data.y)*(atom_data[i].ya-grid_data.y)  +  (atom_data[i].za- grid_data.z)*(atom_data[i].za-grid_data.z));
		if(distance<=mindis)
		{
			minatom=i;
			mindis=distance;
		}
	}
	if(atom_data[minatom].ele==6) return 0;
	else if(atom_data[minatom].ele==1)
	{
		for(j=0;j<= atom_number-1;j++)
		{
			if(minatom!=j)
			{
			
				distance=sqrt(  (atom_data[j].xa-atom_data[minatom].xa)*(atom_data[j].xa-atom_data[minatom].xa) + (atom_data[j].ya-atom_data[minatom].ya)*(atom_data[j].ya-atom_data[minatom].ya)  + (atom_data[j].za-atom_data[minatom].za)*(atom_data[j].za-atom_data[minatom].za) );
				if(distance<1.2)
				{
					if(atom_data[j].ele==7) return 1;
					else if(atom_data[j].ele==8) return 1;
					else return 0;
				}
			
			}
		}
	}
	else if(atom_data[minatom].ele==7) return 1;
	else if(atom_data[minatom].ele==8) return 1;
	else return 0;
		
}










int file_type(char *file_name)
{
	char *mol2=".MOL2";
	char *sdf=".SDF";
	char *pdb=".PDB";
	char *xyz=".XYZ";
	char *p_name;
	char suffix_name[6];
	int i=0;
	p_name=file_name;
	while(*p_name!='\0') p_name++;
	while(*p_name!='.') p_name--;
	while(*p_name!='\0')
	{	
		if(i>=5)
		{	
			printf("\nFile format error\n");
			return -1;
		 } 
		suffix_name[i]=toupper(*p_name);
		p_name++;
		i++;
	} 
	suffix_name[i]='\0';
	if(strcmp(suffix_name, mol2)==0) return 1; 
	else if(strcmp(suffix_name, pdb)==0) return 2; 
	else if(strcmp(suffix_name, sdf)==0) return 3; 
	else if(strcmp(suffix_name, xyz)==0) return 4;
	else
	{
			printf("\nFile format error\n");
			return -1;		
	 } 

}


void usage()
{   
	printf("calmvs [-i<in file name>] [-b < center x, center y, center z, size in the x dimension (Angstroms), size in the y dimension (Angstroms) size in the z dimension (Angstroms)>][Options] [-t<calculate Polar Surface Area>][Options] [-o<out file name>][Options]");
	printf("Try  -H option for more information.");
}
void help()
{   
	printf("calmvs calculates molecule surface volume and compares molecular similarity\n\n\n");
	printf("Usage\n");
	printf("calmvs [-i<in file name>] [-b < center x, center y, center z, size in the x dimension (Angstroms), size in the y dimension (Angstroms) size in the z dimension (Angstroms)>][Options] [-t<calculate Polar Surface Area>][Options] [-o<out file name>][Options]");
	printf("The extension of a file decides the format,Only supports mol2, pdb, sdf formats\n");	
}


int get_atom_info(struct atom_info *atom_info_data)
{
	FILE *fp;
	char fstring[256];
	int atom_number,line_number=0,i,j;
	char *name="ele_data.txt";
	fp=fopen(name,"r");
	
	if(fp!=NULL){
		while(fgets(fstring,255,fp)!=NULL){
			if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    		}
			if(fstring[0]=='#') ;
			else{
				int index; 
				char atom[5],rvdw[10],mw[10];
				sscanf(fstring,"%d %s %s %s",&index,atom,rvdw,mw);
    			char *endptr;
    			float frvdw,fmw;
     			frvdw= strtof(rvdw, &endptr);
				if (*endptr != '\0') { // Check if the entire string has been converted
        			printf("Conversion failed, unconverted parts: %s\n", endptr);
    			} else {
        		;	
    			}
     			fmw= strtof(mw, &endptr);
				if (*endptr != '\0') { // Check if the entire string has been converted
        			printf("Conversion failed, unconverted parts: %s\n", endptr);
    			} else {
        			;
    			}
    			if((0<index)&&(index<199)){
    				atom_info_data[index].falg=1;
    				atom_info_data[index].id=index;
    				atom_info_data[index].mw=fmw;
    				atom_info_data[index].vdwr=frvdw;
    				int ia=0;
    				for(ia=0;ia<=4;ia++) atom_info_data[index].atoms[ia]=atom[ia];
    		
				}else{
					printf("Get atom number data from ele_data.txt failed,  atom number<0 or atom number > 300 \n");
				}
			
			}
		}
					
	}
	else{
		printf ("/n Open file %s failed, please check it\n",name);
		exit(-1);
	}
	return 0;
}


int get_data(struct atoms **atom_data,struct atom_info *atom_info_data,int count, float x,float y,float z,char *ele_s)
{
	
	(*atom_data)[count].xa=x;
	(*atom_data)[count].ya=y;
	(*atom_data)[count].za=z;
	
	if((ele_s[1])=='\0'){
		int temp_i=0,flag=0;
		for(temp_i=0;temp_i<=199;temp_i++){
			if(atom_info_data[temp_i].falg==1){
				if((ele_s[0]==atom_info_data[temp_i].atoms[0])&&(atom_info_data[temp_i].atoms[1]=='\0')){
					(*atom_data)[count].ele=atom_info_data[temp_i].id;
					(*atom_data)[count].rvdw = atom_info_data[temp_i].vdwr;
					flag++;
					break;}else;											
				}else ; 
		}
		if(flag==0){
			printf("\n Error cant get atom type #%s#\n",ele_s);
			return -1;}			
	}else if(ele_s[2]=='\0'){
		int temp_i=0,flag=0 ;
		for(temp_i=0;temp_i<=199;temp_i++){
			if(atom_info_data[temp_i].falg==1){
				if((ele_s[0]==atom_info_data[temp_i].atoms[0])&&(atom_info_data[temp_i].atoms[1]==ele_s[1])){
					(*atom_data)[count].ele=atom_info_data[temp_i].id;
					(*atom_data)[count].rvdw = atom_info_data[temp_i].vdwr;
					flag++;
					break;}												
			}
		}
		if(flag==0){
			printf("\n Error cant get atom type %s\n",ele_s);
			return -1;}									
	}else{
			printf("\nError wrong atom type\n");
			return -1;									
	}	
	
	return 0;
}











int get_data_from_file(char *file_name,struct atom_info *atom_info_data,struct atoms **atom_data,int *p_atom_number, int *xyz_file_number)
{
	FILE *fp;
	int atom_number=0,i=0,j=0,count=0,file_number=0;
	float x,y,z;
	char  fstring[356], ele_s[8],cotain[300],ele[8], ele2[8],atom_arry[5];
	*xyz_file_number=1;
	fp=fopen( file_name ,"r");
	
	if(fp==NULL){
			printf ("/n Open file %s failed\n", file_name );
			return -1;		
	}
	
	if(file_type( file_name )==-1){
		printf ("/n Not supporting file format of file %s \n", file_name );
		return -1;
	}	
	if(file_type( file_name )==1){
		int start=0,end=0,lock_mol2=0,number_mol2,number_mol2_1;
		float charge;
		char new_line[80],mol2_name[15];
		char *atom_start="@<TRIPOS>ATOM",*atom_end="@<TRIPOS>BOND",*linep,*linep1;
		
		if(fp!=NULL){
			while ((fgets(fstring,119,fp))!=NULL){
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if(strncmp(fstring,atom_start,13)==0) start=1;
				if(strncmp(fstring,atom_end,13)==0) end=1;
				if((start==1)&&(end==0)&&(strncmp(fstring,"@",1)!=0)){	
					atom_number++;
				}
			}
		}else{
			printf ("/n Open file %s failed ", file_name );
			fclose(fp);
			return -1;		
		}
		rewind(fp);
		if((atom_number>0)&&(atom_number<=100000)) *atom_data=(struct atoms *)malloc(atom_number*sizeof(struct atoms ));
		else {
			printf ("/n ERROR read the count of atoms failed \n");
			return -1;
		}
		start=0;
		end  =0;
		if(fp!=NULL){
			while ((fgets(fstring,119,fp))!=NULL){
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if(strncmp(fstring,atom_start,13)==0) start=1;
				if(strncmp(fstring,atom_end,13)==0) end=1;
				if((start==1)&&(end==0)&&(strncmp(fstring,"@",1)!=0)){	
              	  linep=fstring;
             	   i=0;
            		while(*linep!='\0'){                  	
                		if(isgraph(*linep)!=0){
							   if(lock_mol2==0){
					   				linep1=linep;
                   					while(isgraph(*linep1)!=0){
                    					cotain[i]=*linep1;
                    					linep1++;
                   						i++;                  				
									}
									cotain[i]=' ';
									i++;						   	
							   }	
                   				lock_mol2=1;
						}
						if(isgraph(*linep)==0) lock_mol2=0;
						linep++;
					}
					cotain[i]='\0';
              		sscanf(cotain,"%d %s %f %f %f %s %d %s %f",&number_mol2,ele,&x,&y,&z,ele2,&number_mol2_1,mol2_name,&charge);
               		i=0;
              		while(ele2[i]!='\0'){
                  	 		ele2[i]=toupper(ele2[i]);
                  	 		if(	ele2[i]=='.'){
                   			ele2[i]='\0';
                   			break;
					  	 } else;
                   		i++;
					}
              		if(count<atom_number){
						get_data(atom_data,atom_info_data,count,x,y,z,ele2);	    
						count++;							
					}else{
						printf("\nError read atom data error\n");
						exit(-1);
					}		
				}
			}
		}else{
			printf ("/n Open file %s failed ", file_name );
			fclose(fp);
			return -1;		
		}
		
		*p_atom_number=atom_number;
	} else if(file_type( file_name )==2)
	{
		char *atom="ATOM",*hetatom="HETATM";	
		if(fp!=NULL){
			while ((fgets(fstring,119,fp))!=NULL){
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if((strncmp(fstring,atom,4)==0)||(strncmp(fstring,hetatom,6)==0)){
					atom_number++;
				}
			}		
		}else{
		printf ("/n Open file %s failed ", file_name );
		fclose(fp);
		return -1;
		}
		rewind(fp);
		if((atom_number>0)&&(atom_number<=100000)) *atom_data=(struct atoms *)malloc(atom_number*sizeof(struct atoms ));
		else {
			printf ("/n ERROR read the count of atoms failed \n");
			return -1;
		}
		if(fp!=NULL){
			while ((fgets(fstring,119,fp))!=NULL){
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if((strncmp(fstring,atom,4)==0)||(strncmp(fstring,hetatom,6)==0)){
					if(fstring[76]==' '){
						ele_s[0]=toupper(fstring[77]);
						ele_s[1]='\0';
					//	printf("%s\n",ele_s);
					}else if(((fstring[76]>='A')&&(fstring[76]<='Z'))||((fstring[76]>='a')&&(fstring[76]<='z'))){
						ele_s[0]=toupper(fstring[76]);
						if(fstring[77]==' ')ele_s[1]='\0';
						else if(((fstring[77]>='A')&&(fstring[77]<='Z'))||((fstring[77]>='a')&&(fstring[77]<='z'))){
							ele_s[1]=toupper(fstring[77]);
							ele_s[2]='\0';							
						}else ele_s[1]='\0';

					//	printf("%s\n",ele_s);						
					}else{
							printf("\nError wrong PDB file data colume 77 and colume 78 no elements data\n");
							return -1;						
					}	
					for(i=0;i<=29;i++) fstring[i]=' ';
					for(i=54;i<=79;i++) fstring[i]=' ';
				//	printf("#%s#\n",ele_s);
					sscanf(fstring,"                                %f  %f  %f                        ",&x,&y,&z);
					if(count<atom_number){
							get_data(atom_data,atom_info_data,count,x,y,z,ele_s);
							count++;
					}
					else{
						printf("\nError read atom data error\n");
						exit(-1);
					}
				}
			}		
		}
		else{
			printf ("/n Open file %s failed \n", file_name );
			fclose(fp);
			return -1;
		}
		
		*p_atom_number=atom_number;
	}else if (file_type( file_name )==3)
	{
		int line_count=0;	
		count=0; 
		if(fp==NULL){
			printf("Open file failed\n");
			return -1;
		}
		while(fgets(fstring,255,fp)!=NULL){  
			if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      			fstring[strlen(fstring)-2]='\n';
      			fstring[strlen(fstring)-1]='\0';
    		} 
			line_count++;	
			if(line_count==4){	
			fstring[3]='\0';
			sscanf(fstring,"%d",&atom_number);
			break;
			} 
		}
		*p_atom_number=atom_number;
		if((atom_number>0)&&(atom_number<=100000)) *atom_data=(struct atoms *)malloc(atom_number*sizeof(struct atoms ));
		else {
			printf ("\n ERROR read the count of atoms failed atom number %d\n",atom_number);
			return -1;
		}
		line_count=1;
		while(fgets(fstring,255,fp)!=NULL){  
		 	if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      			fstring[strlen(fstring)-2]='\n';
      			fstring[strlen(fstring)-1]='\0';
    		}
			if(line_count<=atom_number)	{	
				i=0;
				while(isalpha(fstring[i])==0) i++;
				while(isalpha(fstring[i])) i++;
				fstring[i]='\0';
				i=0; 
				j=0;
				while(isspace(fstring[i])) i++;
				while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
					cotain[j]=fstring[i];
					i++;
					j++;
				}
				cotain[j]=' ';
				j++;
				while(isspace(fstring[i])) i++;
				while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
					cotain[j]=fstring[i];
					i++;
					j++;
				}
				cotain[j]=' ';
				j++;
				while(isspace(fstring[i])) i++;
				while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
					cotain[j]=fstring[i];
					i++;
					j++;
				}
				cotain[j]=' ';
				j++;
				while(isspace(fstring[i])) i++;
				while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
					cotain[j]=fstring[i];
					i++;
					j++;
				}
				cotain[j]='\0';	
				sscanf(cotain,"%f %f %f %s",&x,&y,&z,ele);
				i=0;
           	 while(ele[i]!='\0'){
            	    ele[i]=toupper(ele[i]);
            	    i++;
				}
				if(count<atom_number){
					get_data(atom_data,atom_info_data,count,x,y,z,ele);	    
					count++;							
				}else {
					printf("\nError read atom data error\n");
					exit(-1);
				}		
			}	else break;
			line_count++;		
		}
		
		*p_atom_number=atom_number;
	}	else if (file_type( file_name )==4)
	{
		
		int line_count=0; 
		while(fgets(fstring,255,fp)!=NULL){
			if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      			fstring[strlen(fstring)-2]='\n';
      			fstring[strlen(fstring)-1]='\0';
    		}  	
			line_count++;
			if(line_count==1){
				i=0;
				j=0;
				while(isspace(fstring[i])) i++;
				while(isspace(fstring[i])==0){
					atom_arry[j] =fstring[i];
					i++;
					j++;					
				}
				atom_arry[j]='\0';
				sscanf(atom_arry,"%d",&atom_number);
				
			}
			if(line_count==2) break;
		}
		rewind(fp);
		line_count=0;
		while(fgets(fstring,255,fp)!=NULL){
			if(strlen(fstring)>=2) line_count++;		
		}
		file_number=line_count/(atom_number+2);
		*p_atom_number=atom_number;
		*xyz_file_number=file_number;	
		count=0; 
		rewind(fp);
		if((atom_number>0)&&(atom_number<=100000)) *atom_data=(struct atoms *)malloc(atom_number*file_number*sizeof(struct atoms ));
		else {
			printf ("\n ERROR read the count of atoms failed \n");
			return -1;
		}	
		line_count=0;
		line_count=atom_number+2;
		count=0;
		
		
		while(fgets(fstring,255,fp)!=NULL){
			if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      			fstring[strlen(fstring)-2]='\n';
      			fstring[strlen(fstring)-1]='\0';
    		}
		   
			line_count++;
			if(((line_count-1)%(atom_number+2)!=0)&&((line_count-2)%(atom_number+2)!=0)){
				if(line_count<=(atom_number+2)*(file_number+1)){
					i=0; 
					j=0;				
					if(fstring[0]==' ') while(isspace(fstring[i])) i++;
					while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
						if(isalpha(fstring[i])){
							cotain[j]=fstring[i];
					//	printf("cotain[%d] --[%c]\n",j,cotain[j]);
							i++;       
							j++;							
						}
						else i++;
					}
					cotain[j]=' ';
					j++;
					if(fstring[i]=' ')	while(isspace(fstring[i])) i++;
					while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
						cotain[j]=fstring[i];
						i++;  //x
						j++;
					}
					cotain[j]=' ';
					j++;
					
					while(isspace(fstring[i])) i++;
					
					while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
						cotain[j]=fstring[i];
						i++;   //y
						j++;
					}
					cotain[j]=' ';
					j++;
					while(isspace(fstring[i])) i++;
					while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
						cotain[j]=fstring[i];
						i++;
						j++;
					}
					j++;
					cotain[j]='\0';	
					sscanf(cotain,"%s %f %f %f",ele,&x,&y,&z);
					i=0;
           			while(ele[i]!='\0'){
              			ele[i]=toupper(ele[i]);
               			i++;
					}
					ele[2]='\0';
					if(count<=atom_number*file_number-1){
						get_data(atom_data,atom_info_data,count,x,y,z,ele);	    
						count++;					
					}	
				}else break;			
			}	
		}
		
		
	}else{
		printf("\nError file type error\n");
		fclose(fp);
		exit(-1);		
	}
	fclose(fp);	
	return 	0;		
}

void trace_falef(struct atoms *atom_data_trace,struct atoms *atom_data, int atom_number,int file_number)
{
	int i;
	for(i=0;i<=atom_number-1;i++){
		atom_data[i].xa=atom_data_trace[(file_number-1)*atom_number+i].xa;
		atom_data[i].ya=atom_data_trace[(file_number-1)*atom_number+i].ya;
		atom_data[i].za=atom_data_trace[(file_number-1)*atom_number+i].za;
		atom_data[i].ele=atom_data_trace[(file_number-1)*atom_number+i].ele;
		atom_data[i].rvdw=atom_data_trace[(file_number-1)*atom_number+i].rvdw;
	}
	
return ;
	
	
}


int calculate(struct atoms *atom_data, struct result *resultdata,int if_trace_file,int atom_number,float water_flag,float water_radius,int file_number,float *surface_data,float *psa_data,float *volume_data,int out_data)
{
	float grid_gap;	
	long long int i,j ,k,l, box_x,box_y,box_z,count=0,s_cal_dis2,cal_n,cal_nadd1,cal_nsub1;
	float x,y,z,box_vol;
	int grid_n,line_count=0,*noblank_index;
	struct data_index *data_index;
	int box_a,box_b,box_c,temp;
	struct min_max minmax;
	struct timespec start, end;

    if (clock_gettime(CLOCK_MONOTONIC, &start) == -1) {
        perror("clock_gettime");
        return 1;
    }

	
	
	
	
	
		
	minmax= get_min_max(atom_data,atom_number);
	box_vol=(minmax.xmax-minmax.xmin)*(minmax.ymax-minmax.ymin)*(minmax.zmax-minmax.zmin);
	if(box_vol<100.0) grid_gap=0.40;
	else if ((box_vol>=100.0)&&(box_vol<360.0)) grid_gap=0.44;
	else if((box_vol>=360.0)&&(box_vol<1000.0))	grid_gap=0.42;	
	else grid_gap=0.43;	
	box_x=(int)((minmax.xmax-minmax.xmin+8.0)/grid_gap);
	box_y=(int)((minmax.ymax-minmax.ymin+8.0)/grid_gap);
	box_z=(int)((minmax.zmax-minmax.zmin+8.0)/grid_gap);
	grid_data=(struct grid *)malloc(box_x*box_y*box_z*sizeof(struct grid));
	
  
	data_index=(struct data_index *)malloc(box_x*box_y*box_z*sizeof(struct data_index));




count=0;
x=minmax.xmin-4.0;
for(i=0;i<=(box_x-1);i++)
{	
	y=minmax.ymin-4.0;
	for(j=0;j<=(box_y-1);j++)
	{	
		z=minmax.zmin-4.0;
		for(k=0;k<=(box_z-1);k++)
		{	
			if(count<=(box_x*box_y*box_z))
			{
				grid_data[count].x=x;
				grid_data[count].y=y;
				grid_data[count].z=z;
			
			}
			else break;
			count++;
		

	
			z=z+grid_gap;
		}
		y=y+grid_gap;
	}
	

	x=x+grid_gap;
}



count=(box_x)*(box_y)*(box_z);
for(i=0;i<=count-1;i++) data_index[i].wsp=0;


for(i=0;i<=count-1;i++)
{
grid_data[i].state   =    cal_dis(atom_data,atom_number,water_radius,grid_data[i].x,grid_data[i].y,grid_data[i].z);		
if(grid_data[i].state==-1){
	printf("\nError unsupport atom type.\n");
	return -1;

}	
else {
	if(grid_data[i].state<=2){
		data_index[i].wsp=1;
		data_index[i].x=grid_data[i].x;
		data_index[i].y=grid_data[i].y;
		data_index[i].z=grid_data[i].z;	
	} 
	}
}

int number_of_no_blank=0;
for(i=0;i<=count-1;i++){
	if(data_index[i].wsp==1) number_of_no_blank++;
}

noblank_index=(int *)malloc( number_of_no_blank*sizeof(int));


j=0;
for(i=0;i<=count-1;i++){
	if(data_index[i].wsp==1){
		noblank_index[j]=i;
		j++;
	}
}
for(i=0;i<=count-1;i++)
{

if(grid_data[i].state==100){
grid_data[i].state=0;

}

}



for(i=0;i<=number_of_no_blank-1;i++)
{

grid_data[noblank_index[i]].psa =0;
}




for(i=0;i<=number_of_no_blank-1;i++)
{
//	printf("%d - %d\n",i,cal_adjacent(i, grid_data,(box_x*box_y*box_z),box_x,box_y,box_z));
	if (cal_adjacent(noblank_index[i], grid_data,(box_x*box_y*box_z)-1,box_x,box_y,box_z)==3) 
	{    
		grid_data[noblank_index[i]].state=2;
	
	}
}












for(i=0;i<=number_of_no_blank-1;i++)
{
	
	if (grid_data[noblank_index[i]].state==1) 
	{    
		grid_data[noblank_index[i]].state=2;
	
	}
}
// x
for(i=0;i<=box_z-1;i++)
{
	for(j=0;j<=box_y-1;j++)
	{
		for(k=0;k<=box_x-1;k++)
		{
			cal_n=(k*box_y*box_z)+(j*box_z)+i;
			cal_nadd1=((k+1)*box_y*box_z)+(j*box_z)+i;
			cal_nsub1=((k-1)*box_y*box_z)+(j*box_z)+i;
			if((cal_nsub1>=0)&&(cal_nsub1<=(box_y*box_z*box_x-1)))
			{
				if((cal_nadd1>=0)&&(cal_nadd1<=(box_y*box_z*box_x-1)))
				{
					if((cal_n>=0)&&(cal_n<=(box_y*box_z*box_x-1)))
					{
						if(grid_data[cal_nadd1].state!=1)
						{
							if(grid_data[cal_nsub1].state!=1)
							{
							if(grid_data[cal_nadd1].state != grid_data[cal_nsub1].state ){
								grid_data[cal_n].state=1;
							} 
							
							}
						}
						

					}
				}
			}
			
			
			
			
			
		
		}
	  
	}  
		

}

for(i=0;i<=box_x-1;i++)
{
	for(j=0;j<=box_z-1;j++)
	{
		for(k=0;k<=box_y-1;k++)
		{
			cal_n=(i*box_y*box_z)+(k*box_z)+j;
			
	
			
			
			
			cal_nadd1=(i*box_y*box_z)+((k+1)*box_z)+j;
			cal_nsub1=(i*box_y*box_z)+((k-1)*box_z)+j;

				
			
			
			if((cal_nsub1>=0)&&(cal_nsub1<=(box_y*box_z*box_x-1)))
			{
				if((cal_nadd1>=0)&&(cal_nadd1<=(box_y*box_z*box_x-1)))
				{
					if((cal_n>=0)&&(cal_n<=(box_y*box_z*box_x-1)))
					{
						if(grid_data[cal_nadd1].state!=1)
						{
							if(grid_data[cal_nsub1].state!=1)
							{
								if(grid_data[cal_nadd1].state != grid_data[cal_nsub1].state ) 							
							 		 grid_data[cal_n].state=1;
							}
						}
						

					}
				}
			}
			
			
			
			
			
		
		}
	  
	}  
		

}

//z

for(i=0;i<=box_x-1;i++)
{
	for(j=0;j<=box_y-1;j++)
	{
		for(k=0;k<=box_z-1;k++)
		{
			cal_n=    (i*box_y*box_z)+(j*box_z)+k;
			cal_nadd1=(i*box_y*box_z)+(j*box_z)+k+1;
			cal_nsub1=(i*box_y*box_z)+(j*box_z)+k-1;
	
			if((cal_nsub1>=0)&&(cal_nsub1<=(box_y*box_z*box_x-1)))
			{
				if((cal_nadd1>=0)&&(cal_nadd1<=(box_y*box_z*box_x-1)))
				{
					if((cal_n>=0)&&(cal_n<=(box_y*box_z*box_x-1)))
					{
						if(grid_data[cal_nadd1].state!=1)
						{
							if(grid_data[cal_nsub1].state!=1)
							{
								if(grid_data[cal_nadd1].state != grid_data[cal_nsub1].state ) 
									grid_data[cal_n].state=1;
							}
						}
						

					}
				}
			}
			
			
			
			
			
		
		}
	  
	}  
		

}




count =0;
for(i=0;i<=number_of_no_blank-1;i++)
{
	
	if (grid_data[noblank_index[i]].state==1) 
	{    
		count++;
	
	}
}
float suface;
suface=count*grid_gap;
suface=suface*grid_gap;
suface=1.26100*suface + 13.50600; 
*surface_data=suface;
if(if_trace_file) resultdata[file_number-1].surf=suface;
if(if_trace_file)
(water_flag==0.0) ? printf("%10d %12f",file_number,suface) : printf("%10d %12f",file_number,suface);
else
(water_flag==0.0) ? printf("	suf	%f	",suface) : printf("	SASA_suf	%f	",suface);      

for(i=0;i<=number_of_no_blank-1;i++)
{
	
	if (grid_data[noblank_index[i]].state==1) 
	{    
	grid_data[noblank_index[i]].psa =	cal_psa(grid_data[noblank_index[i]],atom_data,atom_number);
	
	}
}

int psa_count=0;
float psa_suface;


for(i=0;i<=number_of_no_blank-1;i++)
{
	
	if (grid_data[noblank_index[i]].psa==1) 
	{    
	psa_count++;
	
	}
}

psa_suface=((float)psa_count/(float)count)*suface;

*psa_data=psa_suface;
if(if_trace_file)
(water_flag==0) ? printf(" %12f",psa_suface) : printf(" %12f",psa_suface);
else
(water_flag==0) ? printf("	psa	%f	",psa_suface) : printf("	SASA_psa %f	",psa_suface);
if(if_trace_file)
resultdata[file_number-1].psa=psa_suface;





count=0;
for(i=0;i<=number_of_no_blank-1;i++)
{
if((grid_data[noblank_index[i]].state==2)|| (grid_data[noblank_index[i]].state==1))
	{    
		count++;	
	}
}

float volume;
volume=count*grid_gap;
volume=volume*grid_gap;
volume=volume*grid_gap;
*volume_data=volume;
if(if_trace_file)
(water_flag==0) ? printf(" %12f ",volume) : printf(" %12f ",volume);
else
(water_flag==0) ? printf("volume	%f	",volume) : printf("sasa_volume	%f	",volume);
if(if_trace_file)
resultdata[file_number-1].vol=volume;


if (clock_gettime(CLOCK_MONOTONIC, &end) == -1) {
        perror("clock_gettime");
        return 1;
    }
double elapsed_seconds;
 elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

printf(" %f s\n", elapsed_seconds);


if((out_data==1)&&(if_trace_file==0)){

FILE *outsuf_fp;
if((outsuf_fp=fopen("suf_data.XYZ","w"))==NULL) 
{
	printf("error: cant not open file [%s]\n",outsuf_fp);
	return -1;
}

for(i=0;i<=number_of_no_blank-1;i++)
{
	
	if (grid_data[noblank_index[i]].state==1) 
	{    
		fprintf(outsuf_fp,"%.3f %.3f %.3f\n",grid_data[noblank_index[i]].x,grid_data[noblank_index[i]].y,grid_data[noblank_index[i]].z);
	
	}
}



fclose(outsuf_fp);
printf("Saved surface data on suf_data.XYZ\n");
	
}



free(grid_data);

 free(data_index );
 free(noblank_index);

return 0;	
}

int is_obabel_installed() {
    int ret;
#ifdef _WIN32
    // Windows uses the 'where' command to find obabel.exe
    ret = system("where obabel.exe > nul 2>&1");
#else
    // Unix/Linux/MacOS use which command to find obabel
    ret = system("which obabel > /dev/null 2>&1");
#endif
    return (ret == 0);
}

int get_xyz_atom_number(char *name)
{
	FILE *fp;
	char string[300];
	int file_numb;
	fp=fopen(name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",name);
		return -1;
	}
	else{
		if(fgets(string,255,fp)!=NULL){
			if (strlen(string) > 1 && (string[strlen(string)-2]=='\r')){
      			string[strlen(string)-2]='\n';
      			string[strlen(string)-1]='\0';
    		}
			sscanf(string,"%d",&file_numb);
			fclose(fp);
			return file_numb;
		}
		else{
			fclose(fp);
			return -1;
		} 
	}
}

int get_xyz_file_number(char *name)
{
	FILE *fp;
	char string[300];
	int atom_numb=0,linenumber=0;
	atom_numb=get_xyz_atom_number(name);
	if(atom_numb==-1){
		printf("\n\n get %s atom number failed\n",name);
		return -1;
	}
	fp=fopen(name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",name);
		return -1;
	}
	else{
		while(fgets(string,255,fp)!=NULL){
 			linenumber++;
		}
		fclose(fp);
	}
	return linenumber/(atom_numb+2);
}







int cal_confab(char *file_name)
{
	if(is_obabel_installed()!=1){
		printf("The presence of obabel was not detected on this device.\n Please install, recommended version: Open Babel 3.0.0\n");
		exit( -1);	
	}
	int conf=1000,i,conf_number=0,total_conf=0,round=0;
	float rcutoff=0.01,ecutoff=5.0;
	char call_confab[200],string[200],*gen="..generated";
	char number[20],out_name[50];
	
	out_name[0]='\0';
	strcat(out_name,"cos_");
	strcat(out_name,file_name);
	i=0;
	while(out_name[i]!='\0') i++;
	while(out_name[i]!='.') i--;
	out_name[i]='\0';
	strcat(out_name,".xyz");	
	
	remove(out_name);
	remove("confab.log");	
	
	
	
	FILE *log_fp;
	for(rcutoff=0.5;rcutoff>=0.000;rcutoff=rcutoff-0.1) 
	{

			ecutoff=50.0;
		
			round++;
			call_confab[0]='\0';
			strcat(call_confab,"obabel ");
			strcat(call_confab,file_name);
			strcat(call_confab,"  -O  ");
			strcat(call_confab,out_name);
			strcat(call_confab," --confab --verbose ");   //--conf 1000000
			sprintf(number," %.3f ",rcutoff);
			strcat(call_confab," --rcutoff ");
			strcat(call_confab,number);	
			sprintf(number," %.3f ",ecutoff);		
			strcat(call_confab," --ecutoff ");
			strcat(call_confab,number);
			strcat(call_confab,">> confab.log");
			printf("%s\n",call_confab);
			system(call_confab);
			log_fp=fopen("confab.log","r");
			
			if(rcutoff<=0.1)
			{
				
				fclose(log_fp);
				return 0;				
			}
			if(log_fp!=NULL)
			{
				while(fgets(string,199,log_fp)!=NULL)
				{
					if (strlen(string) > 1 && (string[strlen(string)-2]=='\r')){
      					string[strlen(string)-2]='\n';
      					string[strlen(string)-1]='\0';
    				}
					
					if(strncmp(string,"..tot conformations =",21)==0)
					{
						sscanf(string,"..tot conformations = %d",&total_conf);	
					}
					
					if(strncmp(string,gen,11)==0)
					{
						sscanf(string,"..generated %d conformers",&conf_number);
						printf("round %d total conformations %d generation %d conformers ratio %.3f%%\n",round,total_conf,conf_number,((conf_number*100.0)/total_conf));
						if(total_conf<=500){
						if( conf_number==total_conf  )
						{	
							fclose(log_fp);
							return 0;
						}
						else
						{
							fclose(log_fp);
							remove(out_name);
							remove("confab.log");
						}							
							
						}else{
							if( (conf_number>=500)||    (((conf_number*100.0)/total_conf)>=10.0)    )
							{	
								fclose(log_fp);
								return 0;
							}
							else
							{
								fclose(log_fp);
								remove(out_name);
								remove("confab.log");
							}							
							
							
						}

					}

				}
			}
			else{
					printf("\n Open confab.log failed\n") ;
					return -1;
				}
		
		
	}


}


int cal_conformer(char *file_name)
{
	
	//Obabel  9.mol2 -O        9_4.xyz  --conformer  --nconf 200  --writeconformers   --random  
	if(is_obabel_installed()!=1){
		printf("The presence of obabel was not detected on this device.\n Please install, recommended version: Open Babel 3.0.0\n");
		exit( -1);	
	}
	char out_name[50],call_confab[200];
	
	out_name[0]='\0';
	strcat(out_name,"cos_");
	strcat(out_name,file_name);
	int i=0;
	while(out_name[i]!='\0') i++;
	while(out_name[i]!='.') i--;
	out_name[i]='\0';
	strcat(out_name,".xyz");	
	
	remove(out_name);

	call_confab[0]='\0';
	strcat(call_confab,"obabel ");
	strcat(call_confab,file_name);
	strcat(call_confab,"  -O  ");
	strcat(call_confab,out_name);
	strcat(call_confab," --conformer  --nconf 500  --writeconformers   --random  ");
	printf("%s\n",call_confab);
	system(call_confab);
	return 0;		
}

int cal_obabel_mmf(char *file_name)
{
	if(file_type(file_name)!=4){
		printf("\n\n\n    ERROR The input file [%s] for energy minimization should be in xyz format\n\n",file_name);
		exit(-1);
	}

	int atom_number=0,file_number=0;
	int file_count=0;
	char mm_name[100],command[300];
	mm_name[0]='\0';
	
	atom_number=get_xyz_atom_number(file_name);
	file_number=get_xyz_file_number(file_name);



//----------out result-----------
	strcat(mm_name,"mm_");
	strcat(mm_name,file_name);	
	int i=0;
	while(mm_name[i]!='\0') i++;
	while(mm_name[i]!='.') i--;
	mm_name[i]='\0';
	strcat(mm_name,".xyz");
//--------------------------------
	if(file_number==1){
		// Obabel  1.xyz   -O       1_opt.xyz    --minimize  --append "Energy"  --ff MMFF94   --steps 2500
		command[0]='\0';
		strcat(command,"Obabel ");
		strcat(command,file_name);
		strcat(command," -O ");
		strcat(command,mm_name);
		strcat(command," --minimize  --append \"Energy\"  --ff MMFF94   --steps 3000");
		printf("%s",command);
		system(command);			
	}else if(file_number>1){
	//	remove(mm_name);
		FILE *mmresult;
		mmresult=fopen(mm_name, "w");
		if(mmresult==NULL){
			printf("Open file %s fialed\n",mm_name);
			exit(-1);	
		}
//-----------------------METHOD2  trace file-----------------------------------------

	FILE *file;
	char fs[300];
	file=fopen(file_name, "r");
	if(file==NULL){
		printf("Open file %s fialed\n",file_name);
		exit(-1);	
	}
	int cf=0;
	FILE *temp,*mm_temp;
	temp=fopen("temp.xyz","w");
	if(temp==NULL){
		printf("Open file temp.xyz fialed\n");
		fclose(file);
		fclose(mmresult);
		exit(-1);
	}
	
	char fstring[300];	
	while(fgets(fs,255,file)!=NULL){

		
		
		cf++;
		
		
		
		if((cf-((atom_number+2)*file_count))==2){
			fputs("\n",temp);	
		}
		else{
			fputs(fs,temp);
		}
		
		if((cf%(atom_number+2))==0){
			file_count++;
			fclose(temp);
			
			command[0]='\0';
			strcat(command,"Obabel ");
			strcat(command,"temp.xyz");
			strcat(command," -O ");
			strcat(command,"mm_temp.xyz");
			strcat(command," --minimize  --append \"Energy\"  --ff MMFF94   --steps 3000 ");
			printf("Frame %d of %d ",file_count,file_number);
			printf("%s",command);
			system(command);
			FILE *mm_temp;
			mm_temp=fopen("mm_temp.xyz","r");
			if(mm_temp!=NULL){
				while(fgets(fstring,255,mm_temp)!=NULL){
					fputs(fstring, mmresult);
				}
				fclose(mm_temp);
			}
			else{
				printf("Open file mm_temp.xyz fialed\n");
				fclose(mmresult);
				fclose(file);
				exit(-1);
			}
			if(file_count<=(file_number-1))
			temp=fopen("temp.xyz","w");
			if(temp==NULL){
				printf("Open file temp.xyz fialed\n");
				fclose(mmresult);
				fclose(file);
				exit(-1);
			}		
	
		}
		
	}

//-------------------------------------------------------------
		fclose(mmresult);
		remove("mm_temp.xyz");
		remove("temp.xyz");
		
		
	
		
	}
	
	return 0;	
}





int get_mopac_dataf(char *file_name,char *out_name)
{
	char *energy_flag="          FINAL HEAT OF FORMATION =",
		 fstring[200],
		 *zuobiao="                             CARTESIAN COORDINATES",
		 *atom_number="           Empirical Formula: ",
		 atom_arry[5], *mopac_done=" == MOPAC DONE ==";
	FILE *fp,*out_fp;
	int count=0,zuobiao_flag=0,blank_line=0,i=0,atom=0,finished=0;
	float energy_kcal=0.0,energy_kj=0.0,x,y,z;
/*******************************if calculate finished**********************************/	
	fp=fopen(file_name,"r");
		if(fp!=NULL)
	{
		while ((fgets(fstring,199,fp))!=NULL)
			{

				if(strncmp(fstring,mopac_done,17)==0)
				{
					finished=1;
					break;
				}
				
								
			}
	}
	else
	{
		printf ("\n Open file %s failed \n",file_name);
		exit(-1);		
	}
	fclose(fp);	
	if(finished==0) return -1;
	
	
	
	
/*******************************out xyz name**********************************/


	out_fp=fopen(out_name,"w");

/*****************************************************************/	
	
	
	
	
	
	
	
/**********************get atom number ****************************/
	fp=fopen(file_name,"r");
		if(fp!=NULL)
	{
		while ((fgets(fstring,199,fp))!=NULL)
			{
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if(strncmp(fstring,atom_number,29)==0)
				{
					i=0;
					while(fstring[i]!='=')
					{
						fstring[i]=' ';
						i++;
					}
					sscanf(fstring,"                                            = %d atoms",&atom);
					break;
				}
				
								
			}
	}
	else
	{
		printf ("\n Open file %s failed \n",file_name);
		exit(-1);		
	}
	

fprintf(out_fp,"%d\n",atom);
/*************************************************/	
	rewind(fp);
	if(fp!=NULL)
	{
		while ((fgets(fstring,199,fp))!=NULL)
			{
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				
				if(strncmp(fstring,energy_flag,35)==0)
				{
					sscanf(fstring,"          FINAL HEAT OF FORMATION =        %f KCAL/MOL =    %f KJ/MOL",&energy_kcal,&energy_kj);
					fprintf(out_fp,"%s %-15.7f\n","Energy(KJ/MOL)",energy_kj);
				//	printf("%.5f KJ/MOL %.5f KCAL/MOL\n",,energy_kcal);
				}
				if(strncmp(fstring,zuobiao,35)==0)
				{
					zuobiao_flag=1;
				}
				if((zuobiao_flag==1)&&(strncmp(fstring,zuobiao,35)!=0))
				{
					if(strlen(fstring)==1) blank_line++;
					if((strlen(fstring)!=1)&&(blank_line<=1))
					{
						i=0;
						while(isalpha(fstring[i])==0)
						{
							fstring[i]=' ';
							i++;
						}
						sscanf(fstring,"        %s       %f    %f    %f",atom_arry,&x,&y,&z);
						//printf("[a]%s",fstring);
    					fprintf(out_fp,"%3s%15.5f%15.5f%15.5f\n",atom_arry,x,y,z);
  																	
					}

				}								
			}
	}
	else
	{
		printf ("\n Open file %s failed \n",file_name);
		exit(-1);		
	}
	fclose(fp);	
	fclose(out_fp);
	return 0;
}




void print_progress_bar(int progress, int total) {
    int barWidth = 70; 
    float progressPercent = (float)progress / total;
    int pos = barWidth * progressPercent,i;
 
    printf("[");
    for (i = 0; i < barWidth; ++i) {
        if (i < pos) {
            printf("=");
        } else if (i == pos) {
            printf(">");
        } else {
            printf(" ");
        }
    }
    printf("] %d %%\r", (int)(progressPercent * 100));
    fflush(stdout);
}






 
int cal_mopac_mmf(char *file_name,int multiplicity, int charge, float permittivity,char *mopac_path )
{
	if(file_type(file_name)!=4){
		printf("\n\n\n    ERROR The input file for energy minimization should be in xyz format\n\n");
		exit(-1);
	}	

	int atom_number=0,file_number=0;
	
	char mm_name[100],command[300];
	mm_name[0]='\0';
	

	
	atom_number=get_xyz_atom_number(file_name);
	file_number=get_xyz_file_number(file_name);

	
	
	
		
	if(file_number==1){
		
		mm_name[0]='\0';
		strcat(mm_name,"mopac_");
		strcat(mm_name,file_name); 
		int i=0;
		while(mm_name[i]!='\0') i++;
		while(mm_name[i]!='.') i--;
		mm_name[i]='\0';
		strcat(mm_name,".xyz");
		FILE *mopac_temp; 
		mopac_temp=fopen("temp_mopac.mop","w");
		if(mopac_temp==NULL){
			printf("Open file mopac_temp.mop fialed\n");
			exit(-1);			
		}
		else{
			fprintf(mopac_temp,"CHARGE=%d ",charge);
			
			if(multiplicity==1) fprintf(mopac_temp," SINGLET ");
			else if(multiplicity==2) fprintf(mopac_temp," DOUBLET ");
			else if(multiplicity==3) fprintf(mopac_temp," TRIPLET ");
			else if(multiplicity==4) fprintf(mopac_temp," QUARTET ");
			else if(multiplicity==5) fprintf(mopac_temp," QUINTET ");
			else{
				printf("Incorrect multiplicity parameter, please carefully check.\n This method only supports: 1-5.\n");
				fclose(mopac_temp);
				exit(-1);
			}
			
			fprintf(mopac_temp,"PM6-D3H4 precise OPT");
			if((permittivity>1.0)&&(permittivity<200.0))
			fprintf(mopac_temp," eps=%f ",permittivity);
			fprintf(mopac_temp,"\n");
			fprintf(mopac_temp,"molecule\n");
			fprintf(mopac_temp,"All coordinates are Cartesian\n");	
			int i;
			
			FILE *input_filefp; 
			input_filefp=fopen(file_name,"r");
			if(input_filefp==NULL){
				printf("Open file %s fialed\n",file_name);
				exit(-1);			
			}else{
				char fs[300], ele[5];
				int count=0;
				float x,y,z;
				while ((fgets(fs,256,input_filefp))!=NULL){
					if (strlen(fs) > 1 && (fs[strlen(fs)-2]=='\r')){
      					fs[strlen(fs)-2]='\n';
      					fs[strlen(fs)-1]='\0';
    				}
					count++;
					if(count>2){
						sscanf(fs,"%s %f %f %f",ele,&x,&y,&z);
						fprintf(mopac_temp,"%-2s %10.6f  1  %10.6f  1  %10.6f  1  \n",ele,x,y,z);
					}
					
				}
				fclose(input_filefp);
			}			
			fclose(mopac_temp);
			command[0]='\0';
			strcat(command,mopac_path);
			strcat(command," ");
			strcat(command,"temp_mopac.mop");
		//	printf("%s",command);
			system(command);
			if(get_mopac_dataf("temp_mopac.out",mm_name)!=0){
				printf("Failed to run mopac\n");
				exit(-1);	
			}	
		}
		
	}else if(file_number>1){
		
		int i=0,j=0;
		mm_name[0]='\0';
		strcat(mm_name,"mopac_");
		strcat(mm_name,file_name); 
		while(mm_name[i]!='\0') i++;
		while(mm_name[i]!='.') i--;
		mm_name[i]='\0';
		strcat(mm_name,".xyz");		
		FILE *mmresult;
		mmresult=fopen(mm_name,"w");

		if(mmresult==NULL){
			printf("Open file fialed %s\n",mm_name);
			exit(-1);
		}
		
		FILE *mopac_temp; 
		mopac_temp=fopen("temp_mopac.mop","w");
		if(mopac_temp==NULL){
			printf("Open file mopac_temp.mop fialed\n");
			exit(-1);			
		}else{
			fprintf(mopac_temp,"CHARGE=%d ",charge);
			
			if(multiplicity==1) fprintf(mopac_temp," SINGLET ");
			else if(multiplicity==2) fprintf(mopac_temp," DOUBLET ");
			else if(multiplicity==3) fprintf(mopac_temp," TRIPLET ");
			else if(multiplicity==4) fprintf(mopac_temp," QUARTET ");
			else if(multiplicity==5) fprintf(mopac_temp," QUINTET ");
			else{
				printf("Incorrect multiplicity parameter, please carefully check.\n This method only supports: 1-5.\n");
				fclose(mopac_temp);
				exit(-1);
			}
			
			
			
			fprintf(mopac_temp,"PM6-D3H4 precise OPT");
			if((permittivity>1.0)&&(permittivity<200.0))
			fprintf(mopac_temp," eps=%f ",permittivity);
			fprintf(mopac_temp,"\n");
			fprintf(mopac_temp,"molecule\n");
			fprintf(mopac_temp,"All coordinates are Cartesian\n");		
		}
		
		
		FILE *file;
		char fs[300];
		file=fopen(file_name, "r");
		if(file==NULL){
			printf("Open file %s fialed\n",file_name);
			exit(-1);	
		}
		int line=0,file_count=0; 
		float x,y,z;
		char ele[5];
		while(fgets(fs,255,file)!=NULL){
			if (strlen(fs) > 1 && (fs[strlen(fs)-2]=='\r')){
      			fs[strlen(fs)-2]='\n';
      			fs[strlen(fs)-1]='\0';
    		}
			
			line++;
			if((line-((atom_number+2)*file_count))>2){
				sscanf(fs,"%s %f %f %f",ele,&x,&y,&z);
				fprintf(mopac_temp,"%-2s %10.6f  1  %10.6f  1  %10.6f  1  \n",ele,x,y,z);					
			}
			if((line%(atom_number+2))==0){
				file_count++;
				fclose(mopac_temp);
				command[0]='\0';
				
				strcat(command,mopac_path);
				strcat(command," ");
				strcat(command,"temp_mopac.mop");
				printf("%s",command);
				system(command);				
				printf("Running progress:\n");
				print_progress_bar((file_count*100/file_number), 100);
				
				if(get_mopac_dataf("temp_mopac.out","temp_mopac.xyz")!=0){
					printf("Failed to run mopac\n");
					fclose(mmresult);
					exit(-1);	
				}
				FILE *temp_mopac_result;
				temp_mopac_result=fopen("temp_mopac.xyz","r");
				if(temp_mopac_result==NULL){
					printf("Open file fialed temp_mopac.xyz\n");
					fclose(mmresult);
					exit(-1);				
				}else{
					char fs[300];
					while(fgets(fs,255,temp_mopac_result)!=NULL){
						fputs(fs, mmresult);
					}
					fclose(temp_mopac_result);	
					remove("temp_mopac.xyz");
				}				
				
				if(file_count<=(file_number-1)){
					mopac_temp=fopen("temp_mopac.mop","w");
					if(mopac_temp==NULL){
						printf("Open file mopac_temp.mop fialed\n");
						exit(-1);			
					}else{
						fprintf(mopac_temp,"CHARGE=%d ",charge);
						fprintf(mopac_temp,"PM6-D3H4 precise OPT");
						if((permittivity>1.0)&&(permittivity<200.0))
						fprintf(mopac_temp," eps=%f ",permittivity);
						fprintf(mopac_temp,"\n");
						fprintf(mopac_temp,"molecule\n");
						fprintf(mopac_temp,"All coordinates are Cartesian\n");		
					}	
				}

			}

			
		}
		fclose(mmresult);
	
	}else{
		printf("ERROR read file number failed %d\n",file_number);
	}
	
	printf("\n\n %s, file MOPAC energy minimization completed\n",file_name);
	return 0;
}


int get_energy(char *file_name,float bloz)
{
	float blozsel;
	if(bloz==0.0) blozsel=0.5;
	else blozsel=bloz;
	int atom_number,file_number,count=0,fcount=0,linenumber=0;
	FILE *fp;
	atom_number=get_xyz_atom_number(file_name);
	file_number=get_xyz_file_number(file_name);	
	if((atom_number<=1)||(file_number<=1)){
		printf("\n\n read arom number %d or file number %d failed of file %s\n",atom_number,file_number,file_name);
		exit(-1);	
	}
	struct energy{
		int id;
		float energy;
	};
	struct energy  *energy_idx;
	energy_idx=(struct energy *)malloc(file_number*sizeof(struct energy));
	char fs[300];
	fp=fopen(file_name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",file_name);
		return -1;
	}else{
		while(fgets(fs,255,fp)!=NULL){
			if (strlen(fs) > 1 && (fs[strlen(fs)-2]=='\r')){
      			fs[strlen(fs)-2]='\n';
      			fs[strlen(fs)-1]='\0';
    		}
			
			
 			linenumber++;
 			if((linenumber-((atom_number+2)*fcount))==2){
 				float energy;
 				int i=0;
 				while(fs[i]!=' ') {
 					fs[i]=' ';
 					i++;
				 }
 				sscanf(fs,"%f",&energy);
 				energy_idx[fcount].energy=energy;
 				energy_idx[fcount].id=fcount+1;
			 }
			 if((linenumber%(atom_number+2))==0){
			 	fcount++;
			 }
		}	
		fclose(fp);
	}
	int tempid;
	float tempe;
	int i,j;
	for(i=0;i<=file_number-1;i++){
		for(j=i+1;j<=file_number-1;j++){
			if(energy_idx[i].energy>=energy_idx[j].energy){
				tempid=energy_idx[i].id;
				tempe=energy_idx[i].energy;
				energy_idx[i].id=energy_idx[j].id;
				energy_idx[i].energy=energy_idx[j].energy;
				energy_idx[j].id=tempid;
				energy_idx[j].energy=tempe;
			}
		}
	}
	float tempature=298.15,R=0.0019872,dis=0;
	
		for(i=0;i<=file_number-1;i++)
	{
		dis=exp(-(energy_idx[i].energy-energy_idx[0].energy)/(tempature*R))+dis;
	}
	
	printf("\n\nframe  Energy  Distribution(%%)\n");
	
	// choose boltzmann distribution >0.5%
	for(i=0;i<=file_number-1;i++)
	{
		printf("%-5d  %-6.3f  %-6f%%\n",energy_idx[i].id,energy_idx[i].energy,exp(-(energy_idx[i].energy-energy_idx[0].energy)/(tempature*R))/dis*100);
		if((exp(-(energy_idx[i].energy-energy_idx[0].energy)/(tempature*R))/dis*100)<blozsel){
		energy_idx[i].id=-1;	
		}
		
	}
	
	char choosed[100];
	choosed[0]='\0';
	strcat(choosed,file_name);
	i=0;
	while(choosed[i]!='\0') i++;
	while(choosed[i]!='.') i--;
	choosed[i]='\0';
	strcat(choosed,"_choosed.xyz");
	FILE *choose_fp;
	choose_fp=fopen(choosed,"w");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",choosed);
		return -1;
	}
	
	fp=fopen(file_name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",file_name);
		return -1;
	}else{
		
		for(i=0;i<=file_number-1;i++){
			rewind(fp);
			if(energy_idx[i].id>=1){
				linenumber=0;
				while(fgets(fs,255,fp)!=NULL){
 					linenumber++;
 					if((linenumber>((atom_number+2)*(energy_idx[i].id-1)))&&(linenumber<=((atom_number+2)*(energy_idx[i].id)))){
 						fputs(fs,choose_fp);
					 }
				}	
			}
		}
		fclose(fp);
		fclose(choose_fp);
		
	}	
	
	
	
	
	free(energy_idx);
	return 0;
}
 

int get_dist_data(char *file_name,float *distribution, int input_filenumber )
{

	int atom_number,file_number,count=0,fcount=0,linenumber=0;
	FILE *fp;
	atom_number=get_xyz_atom_number(file_name);
	file_number=get_xyz_file_number(file_name);	
	if(input_filenumber!=file_number){
		printf("\n\nThe number of files entered %d does not match the number %d of files in %s file",input_filenumber,file_number,file_name);
		exit(-1);
	}
	if((atom_number<=1)||(file_number<=1)){
		printf("\n\n read arom number %d or file number %d failed of file %s\n",atom_number,file_number,file_name);
		exit(-1);	
	}
	struct energy{
		int id;
		float energy;
	};
	struct energy  *energy_idx;
	energy_idx=(struct energy *)malloc(file_number*sizeof(struct energy));
	char fs[300];
	fp=fopen(file_name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",file_name);
		return -1;
	}else{
		while(fgets(fs,255,fp)!=NULL){
			if (strlen(fs) > 1 && (fs[strlen(fs)-2]=='\r')){
      			fs[strlen(fs)-2]='\n';
      			fs[strlen(fs)-1]='\0';
    		}
			
			
 			linenumber++;
 			if((linenumber-((atom_number+2)*fcount))==2){
 				float energy;
 				int i=0;
 				while(fs[i]!=' ') {
 					fs[i]=' ';
 					i++;
				 }
 				sscanf(fs,"%f",&energy);
 				energy_idx[fcount].energy=energy;
 				energy_idx[fcount].id=fcount+1;
			 }
			 if((linenumber%(atom_number+2))==0){
			 	fcount++;
			 }
		}	
		fclose(fp);
	}
	int tempid;
	float min_energy=energy_idx[0].energy;
	int i;
	for(i=0;i<=file_number-1;i++){
		if(min_energy>=energy_idx[i].energy)
		min_energy=energy_idx[i].energy; 

	}
	float tempature=298.15,R=0.0019872,dis=0;
	
	for(i=0;i<=file_number-1;i++)
	{
		dis=exp(-(energy_idx[i].energy-min_energy)/(tempature*R))+dis;
	}
	
	for(i=0;i<=file_number-1;i++)
	{
		distribution[i]=exp(-(energy_idx[i].energy-min_energy)/(tempature*R))/dis*100;
	}
	free(energy_idx);
	return 0;
}
  


int ask_if_bloz()
{
	int blot;
	char buffer_command[300];
	printf("\nDo we need to perform Boltzmann distribution screening first:enter ""1 (Yes) "" "" 0 (No)""\n" );
	printf("Please enter the selected method:\n");
	fgets(buffer_command, sizeof(buffer_command), stdin);
	sscanf(buffer_command, "%d", &blot);
	return blot;
}
float get_dis()
{
	char buffer_command[300];
	float bloz;
	printf("Enter the lowest proportion of the selected Boltzmann distribution 0-100 (0.5 refers to 0.5%% and default):\n");
	fgets(buffer_command, sizeof(buffer_command), stdin);
	sscanf(buffer_command, "%f", &bloz);
	printf("You selected Boltzmann distribution is£º%f%%\n",bloz);
	return bloz;	
}
void add_prefix(char *newnamef, char *name, char *add)
{
	int i=0;
	newnamef[0]='\0';
	strcat(newnamef,add);
	strcat(newnamef,name);
	while(newnamef[i]!='\0') i++;
	while(newnamef[i]!='.') i--;	
	return ;
}
void add_suffix(char *newnamef, char *name, char *add)
{
	int i=0;
	newnamef[0]='\0';
	strcat(newnamef,name);
	while(newnamef[i]!='\0') i++;
	while(newnamef[i]!='.') i--;
	newnamef[i]='\0';
	strcat(newnamef,add);
	return ;
}



void set_name(int cal_type,char *input_name,char *seted_name)
{
	printf("calculate type is %d\n\n",cal_type);
	char input_file[100],newname[100], newname1[100], newname2[100], newname3[100], newname4[100], newname5[100];
	char name_cal[100];
	int i=0;
	input_file[0]='\0';
	if(cal_type!=1){
		strcat(input_file,input_name);
		i=0;
		while(input_file[i]!='\0'){
			i++;
		}
		while(input_file[i]!='.'){
			i--;
		}	
		input_file[i]='\0';
		strcat(input_file,".xyz");			
	}

	name_cal[0]='\0';
	if(cal_type==1){
		strcat(name_cal,input_name);
	}else if(cal_type==2){  // only conformational search
		add_prefix(name_cal, input_file, "cos_");
	}else if(cal_type==3){    // only MMFF94 bloz
		add_prefix(newname, input_file, "mm_");
		add_suffix(name_cal, newname, "_choosed.xyz");
	}else if(cal_type==4){          // only MMFF94 
		add_prefix(name_cal, input_file, "mm_");
	}else if(cal_type==5){ // only mopac bloz
		add_prefix(newname, input_file, "mopac_");
		add_suffix(name_cal, newname, "_choosed.xyz");		
	}else if(cal_type==6){ // only mopac
		add_prefix(name_cal, input_file, "mopac_");		
	}else if(cal_type==7){ // mmf94 bloz + mopac bloz	
		add_prefix(newname, input_file, "mm_");  // mmff94 mm_a.xyz		
		add_suffix(newname1, newname, "_choosed.xyz"); // bloz mm_a_choose.xyz
		add_prefix(newname2, newname1, "mopac_");  // mopac mopac_mm_a_choose.xyz
		add_suffix(name_cal, newname2, "_choosed.xyz");	// bloz mopac_mm_a_choose_choose.xyz	
	}else if(cal_type==8){
		add_prefix(newname, input_file, "mm_");  // mmff94 mm_a.xyz		
		add_suffix(newname1, newname, "_choosed.xyz"); // bloz mm_a_choose.xyz
		add_prefix(name_cal, newname1, "mopac_");  // mopac mopac_mm_a_choose.xyz		
	}else if(cal_type==9){
		add_prefix(newname, input_file, "mm_");  // mmff94 mm_a.xyz	
		add_prefix(newname1, newname, "mopac_");  // mopac mopac_mm_a.xyz
		add_suffix(name_cal, newname1, "_choosed.xyz"); // bloz mm_a_choose.xyz		
	}else if(cal_type==10){
		add_prefix(newname, input_file, "mm_");  // mmff94 mm_a.xyz	
		add_prefix(name_cal, newname, "mopac_");  // mopac mopac_mm_a.xyz
	}else if(cal_type==11){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_suffix(name_cal, newname1, "_choosed.xyz"); // bloz mm_a_choose.xyz		
	}else if(cal_type==12){
		add_prefix(newname, input_file, "cos_");  // cof cos_a.xyz	
		add_prefix(name_cal, newname, "mm_");  // mmff94 mm_cos_a.xyz
	}else if(cal_type==13){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mopac_");  // mopac mopac_cos_a.xyz
		add_suffix(name_cal, newname1, "_choosed.xyz"); // bloz mopac_cos_a_choose.xyz		
	}else if(cal_type==14){
		add_prefix(newname, input_file, "cos_");  // cof cos_a.xyz	
		add_prefix(name_cal, newname, "mopac_");  // mopac mopac_cos_a.xyz
	}else if(cal_type==15){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_suffix(newname2, newname1, "_choosed.xyz"); // bloz mm_cos_a_choose.xyz	
		add_prefix(newname3, newname2, "mopac_");  // mopac mopac_mm_cos_a_choose.xyz
		add_suffix(name_cal, newname3, "_choosed.xyz"); // bloz mopac_mm_cos_a_choose_choose.xyz		
	}else if(cal_type==16){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_suffix(newname2, newname1, "_choosed.xyz"); // bloz mm_cos_a_choose.xyz	
		add_prefix(name_cal, newname2, "mopac_");  // mopac mopac_mm_cos_a_choose.xyz
	}else if(cal_type==17){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_prefix(newname2, newname1, "mopac_");  // mopac mopac_mm_cos_a.xyz
		add_suffix(name_cal, newname2, "_choosed.xyz"); // bloz mopac_mm_cos_a_choose.xyz		
	}else if(cal_type==18){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_prefix(name_cal, newname1, "mopac_");  // mopac mopac_mm_cos_a.xyz
	}else if(cal_type==19){
		add_suffix(newname, input_file, "_choosed.xyz"); // bloz a_choose.xyz
		add_prefix(newname1, newname, "mopac_");  // mopac mopac_a_choose.xyz
		add_suffix(name_cal, newname1, "_choosed.xyz"); // bloz a_choose_choose.xyz
	}else if(cal_type==20){
		add_prefix(newname, input_file, "mopac_");
		add_suffix(name_cal, newname, "_choosed.xyz");
	}else{
		printf("Error name to calculate\n");
		exit(-1);
	}
	seted_name[0]='\0';
	strcat(seted_name,name_cal);
	
	
	return ;
	
}

void chang_suffix2xyz(char *name)
{
	int i=0;
	while(name[i]!='\0') i++;
	while(name[i]!='.') i--;
	name[i]='\0';
	strcat(name,".xyz");
	return ;
}


// The entire string can be parsed as a number
int is_number(const char *str) {
    char *endptr;
    strtod(str, &endptr);
    return (*endptr == '\0');  
}

