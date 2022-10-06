//================================================================================================//
//
//                       Lattice Boltzmann Method - Functions 3D - Header
//
//================================================================================================//


//	Lê arquivos de inicialização

void read_data ( char*, double*, int*, int*, double*, double*, double* );

//	Lê arquivos de inicialização (modelo térmico com duas distribuições)

void read_data_thermal ( char*, double*, int*, int*, double*, double*, double*, double*, double* );

// Lê arquivo de geometria

int read_geo ( char*, int*, int, int );

// Lê arquivo de geometria (para a determinação da permeabilidade relativa)

int read_geo_kr ( char*, int*, char*, int, int );

// Lê arquivo de geometria acrescentando parede (Zou & He)

int read_geo_ZouHe ( char*, int*, int, int );

// Calcula o produto interno

double prod_int ( double*, double* );

// Arredonda um número

int sum_round ( double  );

// Retorna o momento na direção x

double quant_mov_x ( double*, double* );

// Retorna o momento na direção y

double quant_mov_y ( double*, double* );

// Retorna o momento na direção z

double quant_mov_z ( double*, double* );

// Retorna a vazão mássica em x

double vaz_mass_x ( int*, double*, double*, int, int, int, int );

// Retorna a vazão mássica em y

double vaz_mass_y ( int*, double*, double*, int, int, int, int );

// Retorna a vazão mássica em z

double vaz_mass_z ( int*, double*, double*, int, int, int, int );

// Define os vetores da rede D3Q19

void def_lattice_d3q19 ( double*, double*, double* );

// Define os vetores da rede D3Q15

void def_lattice_d3q15 ( double*, double*, double* );

// Define os vetores da rede D3V59

void def_lattice_d3v59 ( double*, double*, double* );

// Calcula a distribuição de equilíbrio

void dist_eq ( double, double, double, double, double*, double*, double*, double );

// Calcula a distribuição de equilíbrio com Exclusão por volume

void dist_eq_E ( double, double, double, double, double, double, double, double, double, double*, 
				double*, double*, double );

// Calcula a distribuição de equilíbrio  quarta ordem

void dist_eq_thermal ( double, double, double, double, double, double*, double*, double*);

// Calcula densidade e velocidades

void calcula ( double*, double*, double*, double*, double*, double* );

// Calcula o tensor Sum C_alpham C_beta

void calcula_Q ( double*, double* );

// Calcula densidade, velocidades e energia interna

void calcula_thermal ( double*, double*, double*, double*, double*, double*, double* );

// Compara campos de velocidades

double diff_vel_field ( double*, double*, int );

// Calcula a densidade no sítio

double mass ( double* );

// Calcula o momento do sítio

void momentum ( double*, double*, double*, double*, double* );

// Calcula um vetor unitário

void vet_unit ( double, double, double, double*, double*, double* );

// Calcula o operador BGK

void bgk_op ( double*, double*, double, double, double, double, double, double, double, double,
              double*, double*, double );

// Calcula o operador BGK simétrico (modelo TRT)

void bgk_even ( double*, double*, double, double, double, double, double, double, double, double,
                double*, double*, double );

// Calcula o operador BGK simétrico (modelo TRT) com Exclusão por volume

void bgk_even_E ( double*, double*, double, double, double, double, double, double, double, double, 
					double, double, double, double*, double*, double );                

// Calcula o operador BGK antisimétrico (modelo TRT)

void bgk_odd ( double*, double*, double, double, double, double, double, double, double, double,
               double*, double*, double );

// Calcula o operador BGK antisimétrico (modelo TRT) com Exclusão por volume

void bgk_odd_E ( double*, double*, double, double, double, double, double, double, double, double, 
					double, double, double, double*, double*, double );               

// Etapa de colisão usando BGK

void bgk_collision ( double*, double*, double, double, double, double, double*, double );

// Etapa de colisão com regularização - Latt & Chopard 

void regularized_collision ( double*, double*, double*, double, double, double, double );

// Etapa de colisão com regularização - Latt & Chopard + TRT

void regularized_TRT_collision ( double*, double*, double*, double, double );

// Etapa de colisão usando BGK para a rede D3V59

void bgk_collision_thermal ( double*, double*, double, double, double, double, double*, double );

// Etapa de colisão usanto TRT

void trt_collision ( double*, double*, double, double, double, double, double, double*, double );

// Etapa de colisão usanto TRT para distr. theta (advecção/difusão)

void trt_collision_theta ( double*, double*, double, double, double, double,double,double*,double );

// Determina se há variação na posição da interface

double var_inter ( int*, double*, double*, int,	int, int );

// Encontra as direções opostas aos vet. da rede

void i_opposit ( double*, int* );

// Define os sítios para a etapa de propagação

void def_dir_prop ( double*, int*, const int*, int, int, int );

// Define os sítios para a etapa de propagação (com mediadores)

void def_dir_prop_med ( double*, int*, const int*, const bool*, int, int, int );

// Define os sítios para a etapa de propagação otimizada

void def_dir_prop_opt ( double*, int*, int*, int, int, int );

// Etapa de propagação para a distribuição de partículas

void propag_part ( int*, double*, double*, int*, int );

// Etapa de propagação para a distribuição de partículas

void propag_part_opt ( int*, double*, int*, int );

// Etapa de propagação para a distribuição em um único sítio

void propag_part_site ( const int*, double*, double*, int );

// Etapa de propagação para a distribuição g (energia interna) em um único sítio

void propag_g_site ( const int*, double*, double*, const bool*, double, int, double* );

// Etapa de propagação para a distribuição em um único sítio para mediadores

void propag_med_site ( int*, double*, double*, bool*, double, int, double* );

// Determina a posição da interface em y

double pos_interf_y ( int*, double*, double*, int, int, int, int, int, int, int );

// Determina a posição da interface em x

double pos_interf_x ( int*, double*, double*, int, int, int, int, int, int, int );

// Inverte as distribuições apontando para a direção y

void inv_y ( double*, int*, int*, int, int, int, int );

// Inverte as distribuições apontando para a direção z

void inv_z ( double*, int*, int*, int, int, int, int );

// Inverte as distribuições apontando para a direção x

void inv_x ( double*, int*, int*, int, int, int, int );

// Inverte as direções da função distribuição

inline void inv_N ( double* );

// Acrescenta força em um sítio

double force ( double*, double*, double, double, double );

// Grava o meio em formato .vtk

void rec_geometria ( int*, char*, int, int, int );

// Grava os resultados  monofásico

void rec_one_fluid ( int*, double*, double*, double, unsigned int, int, int, int );

// Grava o campo de densidade

void rec_density ( int*, double*, double*, unsigned int, int, int, int, string );

// Grava o campo de velocidades (uma função distribuição)

void rec_velocity ( int*, double*, double*, unsigned int, int, int, int );

// Grava o campo de velocidades (formato binário, uma função distribuição)

void rec_velocity_bin ( int*, double*, double*, unsigned int, int, int, int );

// Grava o campo de velocidades (duas funções distribuição)

void rec_velocity_twofluid ( int*, double*, double*, double*, int, int, int, int );

// Grava os resultados  monofásico térmico (D3Q59)

void rec_one_fluid_thermal ( int*, double*, double*, double, unsigned int, int, int, int );

// Calcula a vorticidade

void vorticity ( int*, double*, double*, int, int, int, int, int, int, double*, double*, double* );

// Grava a vorticidade

void rec_vorticity ( int*, double*, double*, unsigned int, int, int, int );

// Grava arquivo de recuperação

void rec_recovery ( double*, double*, int, int, int );

// Grava arquivo de recuperação (dois fluidos)

void rec_recovery_two_fluids ( double*, double*, double*, int, int );

// Condição de contorno de derivada nula da velocidae direção x

void bond_eqdvnull_x ( int*, double*, double*, int, int, int, int, double*, double );

// Condição de contorno de derivada nula da velocidade mantendo a densidade x

void bond_eqdvnull_rholocal_x ( int*, double*, double*, int, int, int, int, double*, double );

// Condição de contorno de derivada nula da velocidade mantendo a densidade y

void bond_eqdvnull_rholocal_y ( int*, double*, double*, int, int, int, int, double*, double );

// Condição de contorno de derivada nula da velocidade com atenuação, rho local x

void bond_eqdvnull_rholocalveldim_x ( int*, double*, double*, int, double, int, int, int, double*,
                                      double );

// Condição de contorno de derivada nula da velocidade com atenuação, rho local x

void bond_eqdvnull_rholocalveldim_y ( int*, double*, double*, int, double, int, int, int, double*,
                                      double );

// Condição de contorno de derivada nula (Neumann)  direção x (Sem equilíbrio)

void bond_dvnull_rho_x ( int*, double*, double*, int, double, int, int, int, double*, double );

// Condição de contorno de derivada nula (Neumann)  direção x (Sem equilíbrio - d3q19)

void bond_dvnull_x_d3q19 ( int*, double*, double*, int, int, int, int, double*, double );

// Condição de contorno de convectiva (cbc) direção x

void bond_cbc_x ( int*, double*, double*, double*, int, double, int, int, int, double*, double );

// Condição de contorno de convectiva (cbc) direção x para dois fluidos

void bond_cbc_2fase_x ( int*, double*, double*, double*, double*, int, double, int, int, int, 
						double*, double );

// Condição de contorno de convectiva (cbc) direção x para mediadores

void bond_cbc_med_x ( int*, double*, double*, double*, double*, double*, int, double, int, int, 
						int, double*, double );						

// Condição de contorno Zou & He para rede D3Q19 - direção x

void bond_ZouHeD3Q19_rhovel_x ( int* , double*, int, double, double, double, double, int, int, int,
                                double*, double );

// Condição de contorno Zou & He para rede D3Q19, impondo somente rho - direção x

void bond_ZouHeD3Q19_rho_x ( int*, double*, double*, int, double, int, int, int, double*, double );

// Condição de contorno Zou & He para rede D3Q19, impondo somente vx - direção x

void bond_ZouHeD3Q19_vx_x ( int*, double*, double*, int, double, int, int, int, double*, double );

// Condição de contorno de equilíbrio com realimentação (rho & vel.) direção x

void bond_eqretro_rhovel_x ( int*, double*, double*, int, double, double, double, double, int, int,
                             int, double*, double );

// Condição de contorno de equilíbrio com realimentação (rho & vel.) direção y

void bond_eqretro_rhovel_y ( int*, double*, double*, int, double, double, double, double, int, int,
                             int, double*, double );
// Condição de reflexação expecular no plano yz

void bond_mirror_x ( double*, int*, double*, int, int, int, int );

// Condição de reflexação expecular no plano xz

void bond_mirror_y ( double*, int*, double*, int, int, int, int );

// Condição de reflexação expecular no plano xy

void bond_mirror_z ( double*, int*, double*, int, int, int, int );

// Condição de contorno de derivada nula (Neumann)  direção y

void bond_eqdvnull_rho_y ( int*, double*, double*, int, double, int ,int ,int, double*, double );

// Condição de contorno de derivada nula (Neumann)  direção x

void bond_eqdvnull_vx_x ( int*, double*, double*, int, double, int ,int ,int, double*, double );

// Condição de contorno de derivada nula (Neumann)  direção y

void bond_eqdvnull_vy_y ( int*, double*, double*, int, double, int ,int ,int, double*, double );


// Condição de contorno de equilíbrio  direção x

void bond_eq_rhovel_x ( int*, double*, double*, int, double, double, double, double, int, int, int,
                        double*, double );

// Condição de contorno de equilíbrio  direção y

void bond_eq_rhovel_y ( int*, double*, double*, int, double, double, double, double, int, int, int,
                        double*, double );

// Condição de contorno de derivada nula (Neumann)  direção z

void bond_eqdvnull_rho_z ( int*, double*,  double*, int, double, int, int, int, double*, double );

// Condição de contorno de equilíbrio  direção z

void bond_eq_rhovel_z ( int*, double*, double*, int, double, double, double, double, int ,int ,int,
                        double*, double );

// Condição de contorno parabólica  direção x  parede em y

void bond_eqparabolly_x ( int*, double*, double*, int, double, double, int, int, int, int,
                          double*, double );

// Condição de contorno parabólica  direção x  parede em y (Zou & He)

void bond_parabolly_x ( int*, double*, double*, int, double, int, int, int, int,
                          double*, double ); 
           
//	Calcula a permeabilidade intrínseca

double intrinsic_permeability ( double, double, double, double, double, int, int, int, int  );

// Lê arquivo de inicialização (para dois fluidos)

void read_data_twofluid ( char*, double*, int*, int*, double*, double*, double*, double*, double*,
                          double*, double*, double*, double*, double*, double*, string* );

// Calcula um gradiente de um campo escalar

void gradient ( double*, double*, double*, double*, double*, int, int, int, int, int, int, double*,
                double );
                
// Impõe condição de contorno para cálculo da permeabilidade relativa                

void bound_kr ( int*, double*, double*, double*, int, int, int, double*, double*, double*, int, int, 
				int, double*, double );             
                
// Calcula um gradiente de um campo escalar refletindo em um plano em x ou y

void gradient_with_mirror ( double*, double*, int, int, int, int, int, int, int, int, double*, 
							double );

// Calcula um gradiente de um campo escalar refletindo nos sólidos

void gradient_mirror ( double*, double*, int*, double*, double*, double*, int, int, int, int, int, 
						int, double*, double );    

// Calcula um gradiente de um campo escalar usando o ponto no lugar dos sólidos

void gradient_point ( double*, double*, int*, double*, double*, double*, int, int, int, int, int, 
						int, double*, double );    						            

// Etapa de recoloração (Latva-Koko)

void recolloring ( double*, double*, double*, double*, double, double, double, double );

// Retorna a magnitude do campo de velocidades 

void return_vel_field ( double*, double*, double*, double*, int );

// Imposição de tensão interfacial (Spencer, Halliday & Care)

void imp_interf_tension ( double*, double*, double, double, double, double, double*, double );

// Calcula a "densidade efetiva"

double effetive_mass ( double*, double );

// Etapa de colisão usando o modelo Shan & Chen

void sc_collision ( double*, double*, double, double, double, double, double, double, double,
                    double, double*, int, int, int, int, int, int );

// Etapa de colisão usando o modelo Santos, Facin & Philippi

void sfp_collision ( double*, double*, double*, double, double, double, double, double, double,
                     double, double, double, double, double, double, double, double*, double );

// Etapa de colisão usando o modelo Santos, Facin & Philippi com dois tempos de relaxação

void sfp_collision_TRT ( double*, double*, double*, double, double, double, double, double, double,
						 double, double, double, double, double, double, double, double, double, 
						 double*, double );
						 
// Etapa de colisão usando o modelo Santos, Facin & Philippi (TRT) e diferentes densidades

void sfp_collision_TRT_two_factors ( double*, double*, double*, double, double, double, double, 
							double, double, double, double, double, double, double, double, double,
							double, double, double, double*, double );

// Etapa de colisão usando o modelo Santos, Facin & Philippi (TRT) e exclusão por volume

void sfp_collision_TRT_exc ( double*, double*, double*, double*, double, double, double, double, 
							double, double, double, double, double, double, double, double, double, 
							double, double, double*, double, int, int, int, int, int, int );

// Etapa de colisão usando o modelo Santos, Facin & Philippi (TRT) e exclusão por volume

void sfp_collision_TRT_exc_point ( double*, double*, double*, double*, int*, double, double, double, 
							double, double, double, double, double, double, double, double, double,
							double, double, double, double*, double, int, int, int, int, int, int );	

// Etapa de colisão usando o modelo Santos, Facin & Philippi (TRT) e exclusão por volume

void sfp_collision_TRT_exc_point_recoll ( double*, double*, double*, double*, int*, double, double, 
							double, double,
							double, double, double, double, double, double, double, double, double,
							double, double, double, double*, double, int, int, int, int, int, int );													 

// Impõe uma condição inicial aleatória

void ini_cond_rand ( double*, int*, double*, double*, double, double, double, int, int, int,
                     double*, double );

// Impõe uma condição inicial via arquivos

void ini_cond_file ( double*, int*, double*, double*, double, double, string, double, double,
                     double, double*, double );

// Retorna o número de pontos de um dos fluidos

int n_ptos ( double*, double*, int );

// Retorna o número de pontos de um dos fluido (parte da geometria)

int n_ptos_part ( double*, double*, int*, int, int, int, int, int, int, int, int, int );

// Etapa de propagação para a distribuição de mediadores

double propag_med ( int*, double*, double*, int*, char*, double, int );

// Grava os resultados - dois fluidos

void rec_two_fluid ( int*, double*, double*, double*, unsigned int, int, int, int );

// Grava os resultados - dois fluidos com exclusão por volume

void rec_two_fluid_exc ( int*, double*, double*, double*, double*, unsigned int, int, int, int,
                         double*, double);

// Grava os resultados - dois fluidos (Shan & Chen)

void rec_two_fluid_SC ( int*, double*, double*, double*, double*, double*, unsigned int, int,
                        int, int );

// Grava a interface entre dois fluidos

void rec_interface ( int*, double*, double*, unsigned int, int, int, int );

// Etapa de colisão immiscível usando o modelo SHC (Phil. Trans.)

void shc_coll2011 ( double*, double*, double*, double*, double, double, double, double, double,
                         double, double, double, double, double, double, double, double*, double );
                         
// Colisão immiscível, modelo SHC (Phys. Rev.), One Relax. Times                
                         
void shc_coll2010 ( double*, double*, double*, double, double, double, double, double, double, 
					double, double, double, double, double, double, double*, double );                         
    

// Colisão immiscível, modelo SHC (Phys. Rev.), Two Relax. Times

void shc_coll2010_TRT ( double*, double*, double*, double, double, double, double, double, double,
                             double, double, double, double, double,
                             double, double, double*, double );

// Calcula o termo de exclusão por volume

void omega_exc ( double, double, double, double, double, double, double, double*, double*, double*,
                 double );
                 
// Calcula o termo de exclusão por volume MODEL A

void omega_exc_A ( double, double, double*, double*, double*, double );

// Colisão immiscível,modelo SHC (Phys. Rev.), com exclusão

void shc_coll2010_TRT_exc ( double*, double*, double*, double*, double, double, double, double,
                            double, double, double, double, double, double, double, double, double,
                            int, int, int, int, int, int, double*, double );

// Colisão immiscível,modelo SHC (Phys. Rev.), com exclusão MODELO E (Eduardo)

void shc_coll2010_TRT_exc_E ( double*, double*, double*, double*, double, double, double, double,
                            double, double, double, double, double, double, double, double, double,
                            int, int, int, int, int, int, double*, double );                            
                            
// Colisão immiscível,modelo SHC (Phys. Rev.), com exclusão MODEL A

void shc_coll2010_TRT_exc_A ( double*, double*, double*, double*, int*, double, double, double, 
									double, double, double, double, double, double, double, double, 
									double, double, int, int, int, int, int, int, double*, double );

// Colisão immiscível,modelo SHC (Phys. Rev.), com exclusão

void shc_coll2010_TRT_exc_mirror ( double*, double*, double*, double*, int*, double, double, double, 
									double, double, double, double, double, double, double, double, 
									double, double, int, int, int, int, int, int, double*, double );                            

// Colisão immiscível,modelo SHC (Phys. Rev.), com exclusão

void shc_coll2010_TRT_exc_point ( double*, double*, double*, double*, int*, double, double, double, 
									double, double, double, double, double, double, double, double, 
									double, double, int, int, int, int, int, int, double*, double ); 
									
// Calcula a tensão interfacial

double calc_inter_tension ( double*, double*, int*, int, int, int );

//================================================================================================//







//================================================================================================//
//
//                      Lattice Boltzmann Method - Functions 3D
//
//================================================================================================//



//===================== Lê o arquivo de inicialização ============================================//
//
//      Input:
//      Output: nome do arquivo de geometria, tam. do pixel, passos, files, tempo de relaxação,
//              densidade inicial
//
//================================================================================================//

void read_data ( char *nome_geo, double *ftesc, int *npassos, int *numarq, double *tau,
                 double *visc, double *ro_ini )

{
    //--------------------------------------------------------------------------------------------//

    char nome_in[] = "data_in.txt";

    ifstream f_in( nome_in );

    //--------------------------------------------------------------------------------------------//

    char nome_out[] = "dat_out.txt";

    cout << "\nNome do arquivo de dados (saida): " << nome_out << endl;

    //--------------------------------------------------------------------------------------------//

    ofstream fdat( nome_out );

    fdat << "Nome do arquivo de dados (saida): " << nome_out << endl;

    //--------------------------------------------------------------------------------------------//

    f_in >> nome_geo;

    cout << "\nNome do arquivo de geometria: " << nome_geo << endl;

    fdat << "Nome do arquivo de geometria: " << nome_geo << endl;

    //--------------------------------------------------------------------------------------------//

    char st_ft[25];
    f_in >> st_ft;

    f_in >> *ftesc;  // Le a dimensao do dimensão do pixel ( m )

    cout << "\n\nDimensao do pixel = " << *ftesc << endl;

    fdat << "\nDimensao do pixel = " << *ftesc << endl;

    //--------------------------------------------------------------------------------------------//

    char st_steps[25];
    f_in >> st_steps;

    f_in >> *npassos;

    cout << "\nNumero de passos: " << *npassos << endl;

    fdat << "Numero de passos: " << *npassos << endl;

    //--------------------------------------------------------------------------------------------//

    char st_files[25];
    f_in >> st_files;

    f_in >> *numarq;

    cout << "Numero de arquivos: " << *numarq << endl;

    fdat << "Numero de arquivos: " << *numarq << endl;

    //--------------------------------------------------------------------------------------------//

    char st_tau[25];
    f_in >> st_tau;

    f_in >> *tau;

    cout << "Tempo de relaxacao: " << *tau << endl;

    fdat << "Tempo de relaxacao: " << *tau << endl;

    *visc = 1.0 / 3.0 * ( *tau - 0.5 );

    fdat << "viscosidade = " << *visc << endl;

    //--------------------------------------------------------------------------------------------//

    char st_ro[25];
    f_in >> st_ro;

    f_in >> *ro_ini;

    cout << "Densidade: " << *ro_ini << endl;

    fdat << "Densidade: " << *ro_ini << endl;

    //--------------------------------------------------------------------------------------------//

    f_in.close();

    fdat.close();
}

//================================================================================================//




//===================== Lê o arquivo de inicialização (inclui difusão/advecção de calor) =========//
//
//      Input:
//      Output: nome do arquivo de geometria, tam. do pixel, passos, files, tempo de rel. do fluido,
//       viscosidade, tempo de rel.o da temp., coef. de dif de calor, densidade inicial
//
//================================================================================================//

void read_data_thermal ( char *nome_geo, double *ftesc, int *npassos, int *numarq, double *tau_f,
                 double *visc, double *tau_q, double *coef_diff, double *ro_ini )

{
    //--------------------------------------------------------------------------------------------//

    char nome_in[] = "data_in.txt";

    ifstream f_in( nome_in );

    //--------------------------------------------------------------------------------------------//

    char nome_out[] = "dat_out.txt";

    cout << "\nNome do arquivo de dados (saida): " << nome_out << endl;

    //--------------------------------------------------------------------------------------------//

    ofstream fdat( nome_out );

    fdat << "Nome do arquivo de dados (saida): " << nome_out << endl;

    //--------------------------------------------------------------------------------------------//

    f_in >> nome_geo;

    cout << "\nNome do arquivo de geometria: " << nome_geo << endl;

    fdat << "Nome do arquivo de geometria: " << nome_geo << endl;

    //--------------------------------------------------------------------------------------------//

    char st_ft[25];
    f_in >> st_ft;

    f_in >> *ftesc;  // Le a dimensao do dimensão do pixel ( m )

    cout << "\n\nDimensao do pixel = " << *ftesc << endl;

    fdat << "\nDimensao do pixel = " << *ftesc << endl;

    //--------------------------------------------------------------------------------------------//

    char st_steps[25];
    f_in >> st_steps;

    f_in >> *npassos;

    cout << "\nNumero de passos: " << *npassos << endl;

    fdat << "Numero de passos: " << *npassos << endl;

    //--------------------------------------------------------------------------------------------//

    char st_files[25];
    f_in >> st_files;

    f_in >> *numarq;

    cout << "Numero de arquivos: " << *numarq << endl;

    fdat << "Numero de arquivos: " << *numarq << endl;

    //--------------------------------------------------------------------------------------------//

    char st_tau_f[25];
    f_in >> st_tau_f;

    f_in >> *tau_f;

    cout << "Tempo de relaxacao do fluido: " << *tau_f << endl;

    fdat << "Tempo de relaxacao do fluido: " << *tau_f << endl;

    *visc = 1.0 / 3.0 * ( *tau_f - 0.5 );

    fdat << "viscosidade = " << *visc << endl;
    
    //--------------------------------------------------------------------------------------------//

    char st_tau_q[25];
    f_in >> st_tau_q;

    f_in >> *tau_q;

    cout << "Tempo de relaxacao para temperatura: " << *tau_q << endl;

    fdat << "Tempo de relaxacao para temperatura: " << *tau_q << endl;

    *coef_diff = 1.0 / 3.0 * ( *tau_q - 0.5 );

    fdat << "Coeficiente de difusao de calor = " << *coef_diff << endl;

    //--------------------------------------------------------------------------------------------//

    char st_ro[25];
    f_in >> st_ro;

    f_in >> *ro_ini;

    cout << "Densidade: " << *ro_ini << endl;

    fdat << "Densidade: " << *ro_ini << endl;

    //--------------------------------------------------------------------------------------------//

    f_in.close();

    fdat.close();
}

//================================================================================================//




//===================== Lê o arquivo de geometria ================================================//
//
//      Input: nome do arquivo de geometria, ponteiro para o meio
//      Output: número de pontos lidos
//
//================================================================================================//

int read_geo ( char *nome_geo, int *meio, int pts_in, int pts_out )
{
    int nx, ny, nz;

    ifstream fmatriz( nome_geo );
    
	string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;

    cout << "\nDimensao x = " << nx << endl;
    
    cout << "Dimensao y = "   << ny << endl;
    
    cout << "Dimensao z = "   << nz << endl;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    //--------------- Acrescenta layers ----------------------------------------------------------//

    nx = nx + pts_in + pts_out;

    //--------------------------------------------------------------------------------------------//

    int poros = 0;
  
    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;

                if ( x >= pts_in && x < nx - pts_out )
                {
                    fmatriz >> meio[pos];

                    if ( meio[pos] > 0 )
                    {
                        poros++;

                        meio[pos] = poros;
                    }
                }
                else
                {
                    poros++;

                    meio[pos] = poros;
                }
            }
        }
    }

    fmatriz.close();

    return poros;
}

//================================================================================================//




//===================== Lê o arquivo de geometria (para a determinação de Kr) ====================//
//
//      Input: nome do arquivo de geometria, ponteiro para o meio
//      Output: número de pontos lidos
//
//================================================================================================//

int read_geo_kr ( char *nome_geo, int *meio, char* interface, int pts_in, int pts_out )
{
    int nx, ny, nz;

    ifstream fmatriz( nome_geo );
    
	string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;

    cout << "\nDimensao x = " << nx << endl;
    
    cout << "Dimensao y = "   << ny << endl;
    
    cout << "Dimensao z = "   << nz << endl;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    //--------------- Acrescenta layers ----------------------------------------------------------//

    nx = nx + pts_in + pts_out;

    //--------------------------------------------------------------------------------------------//

    int poros = 0;
  
    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;

                if ( x >= pts_in && x < nx - pts_out )
                {
                    fmatriz >> meio[pos];

                    if ( meio[pos] != 0 )
                    {
                        if ( meio[pos] < 0 ) interface[pos] = 1;
                        
						poros++;

						meio[pos] = poros;      
                    }
                }
                else
                {
                    poros++;

                    meio[pos] = poros;
                }
            }
        }
    }

    fmatriz.close();

    return poros;
}

//================================================================================================//




//===================== Lê o arquivo de geometria ================================================//
//
//      Input: nome do arquivo de geometria, ponteiro para o meio
//      Output: número de pontos lidos
//
//================================================================================================//

int read_geo_ZouHe ( char *nome_geo, int *meio, int pts_in, int pts_out )
{
    int nx, ny, nz;

    ifstream fmatriz( nome_geo );
    
	string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( fmatriz, dump );
	
	getline( fmatriz, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;

    cout << "\nDimensao x = " << nx << endl;
    cout << "Dimensao y = "   << ny << endl;
    cout << "Dimensao z = "   << nz << endl;
    
    for ( int i = 0; i < 5; i++ ) getline( fmatriz, dump );

	dados.clear();

    //--------------- Acrescenta layers ----------------------------------------------------------//

    nx = nx + pts_in + pts_out;

    //--------------------------------------------------------------------------------------------//

    int poros = 0;
  
    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int pos = x + y * nx + z * ny * nx;

                if ( x >= pts_in && x < nx - pts_out )
                {
                    fmatriz >> meio[pos];

                    if ( meio[pos] )
                    {
                        poros++;

                        meio[pos] = poros;
                    }
                }
                
                else if ( x == 0 )
                {
                    meio[pos] = 0;
                }
                
                else if ( x == nx - 1 )
                {
                    meio[pos] = 0;
                }
                
                else
                {
                    poros++;

                    meio[pos] = poros;
                }
            }
        }
    }

    fmatriz.close();

    return poros;
}

//================================================================================================//




//===================== Calcula o produto interno ================================================//
//
//      Input: two vetors
//      Output: internal product
//
//================================================================================================//

double  prod_int ( double  vet_1[dim], double  vet_2[dim] )
{

    double  result = 0.0;

    for ( int i = 0; i < dim; i++ )
    {
        result = result + vet_1[i] * vet_2[i];
    }

    return result;
}

//================================================================================================//




//===================== Arredonda um número  =====================================================//
//
//      Input: double
//      Retorna: int
//
//================================================================================================//

int  sum_round ( double  num )
{
    int num_int;

    if ( num < 0 ) num_int = ( int ) ( num - 0.5 );
    else num_int = ( int ) ( num + 0.5 );

    return num_int;
}

//================================================================================================//




//===================== Retorna o momento na direção x ===========================================//
//
//      Input: distribution function *f, lattice vectors
//      Output: momentum in the x direction
//
//================================================================================================//

double quant_mov_x ( double *f, double *ini_c )
{
    double *c;

    double mx = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        mx = mx + c[0] * f[i];
    }

    return mx;

}

//================================================================================================//




//===================== Retorna o momento na direção y ===========================================//
//
//      Input: distribution function *f, lattice vectors
//      Output: momentum in the x direction
//
//================================================================================================//

double quant_mov_y ( double *f, double *ini_c )
{
    double *c;

    double my = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        my = my + c[1] * f[i];
    }

    return my;

}

//================================================================================================//




//===================== Retorna o momento na direção z ===========================================//
//
//      Input: distribution function *f, lattice vectors
//      Output: momentum in the z direction
//
//================================================================================================//

double quant_mov_z ( double *f, double *ini_c )
{
    double *c;

    double mz = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        mz = mz + c[2] * f[i];
    }

    return mz;

}

//================================================================================================//




//===================== Retorna a vazão mássica em x =============================================//
//
//      Input: geometry, distribution function *f, lattice vectors
//      Output: mass flow
//
//================================================================================================//

double vaz_mass_x ( int *ini_meio, double *ini_f, double *ini_c, int pos_x, int nx, int ny, int nz )
{
    int *meio;

    double sum_mx = 0.0;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            meio = ini_meio + pos_x + y * nx + z * ny * nx;

            if ( *meio )
            {
                double *f;

                f = ini_f + ( *meio - 1 ) * nvel;

                sum_mx = sum_mx + quant_mov_x ( f, ini_c );
            }
        }
    }

    return sum_mx;
}

//================================================================================================//





//===================== Retorna a vazão mássica em y =============================================//
//
//      Input: geometry, distribution function *f, lattice vectors
//      Output: mass flow
//
//================================================================================================//

double vaz_mass_y ( int *ini_meio, double *ini_f, double *ini_c, int pos_y, int nx, int ny, int nz )
{
    int *meio;

    double sum_my = 0.0;

    for ( int x = 0; x < nx; x++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            meio = ini_meio + x + pos_y * nx + z * ny * nx;

            if ( *meio )
            {
                double *f;

                f = ini_f + ( *meio - 1 ) * nvel;

                sum_my = sum_my + quant_mov_y ( f, ini_c );
            }
        }
    }

    return sum_my;
}

//================================================================================================//




//===================== Retorna a vazão mássica em z =============================================//
//
//      Input: geometry, distribution function *f, lattice vectors
//      Output: mass flow
//
//================================================================================================//

double vaz_mass_z ( int *ini_meio, double *ini_f, double *ini_c, int pos_z, int nx, int ny, int nz )
{
    int *meio;

    double sum_mz = 0.0;

    for ( int x = 0; x < nx; x++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            meio = ini_meio + x + y * nx + pos_z * ny * nx;

            if ( *meio )
            {
                double *f;

                f = ini_f + ( *meio - 1 ) * nvel;

                sum_mz = sum_mz + quant_mov_z ( f, ini_c );
            }
        }
    }

    return sum_mz;
}

//================================================================================================//




//===================== Define os vetores da rede D3Q19 ==========================================//
//
//      Input: pointer to the vectors
//      Output:
//
//================================================================================================//

void def_lattice_d3q19 ( double *ini_c, double *W, double *one_over_c_s2 )
{
    double  *c;

    //------------- |c| = 0 ----------------------//

    c = ini_c;
    c[0] = 0;
    c[1] =  0;
    c[2] = 0;

    //------------- |c| = 1 ----------------------//

    c = ini_c + 1 * dim;
    c[0] =  1;
    c[1] =  0;
    c[2] =  0;

    c = ini_c + 2 * dim;
    c[0] = -1;
    c[1] =  0;
    c[2] =  0;

    c = ini_c + 3 * dim;
    c[0] =  0;
    c[1] =  1;
    c[2] =  0;

    c = ini_c + 4 * dim;
    c[0] =  0;
    c[1] = -1;
    c[2] =  0;

    c = ini_c + 5 * dim;
    c[0] =  0;
    c[1] =  0;
    c[2] =  1;

    c = ini_c + 6 * dim;
    c[0] =  0;
    c[1] =  0;
    c[2] = -1;

    //------------ |c| = sqrt(2) -----------------//

    c = ini_c + 7 * dim;
    c[0] =  1;
    c[1] =  1;
    c[2] =  0;

    c = ini_c + 8 * dim;
    c[0] = -1;
    c[1] = -1;
    c[2] =  0;

    c = ini_c + 9 * dim;
    c[0] =  1;
    c[1] = -1;
    c[2] =  0;

    c = ini_c + 10 * dim;
    c[0] = -1;
    c[1] =  1;
    c[2] =  0;

    c = ini_c + 11 * dim;
    c[0] =  1;
    c[1] =  0;
    c[2] =  1;

    c = ini_c + 12 * dim;
    c[0] = -1;
    c[1] =  0;
    c[2] = -1;

    c = ini_c + 13 * dim;
    c[0] =  1;
    c[1] =  0;
    c[2] = -1;

    c = ini_c + 14 * dim;
    c[0] = -1;
    c[1] =  0;
    c[2] =  1;

    c = ini_c + 15 * dim;
    c[0] =  0;
    c[1] = -1;
    c[2] = -1;

    c = ini_c + 16 * dim;
    c[0] =  0;
    c[1] =  1;
    c[2] =  1;

    c = ini_c + 17 * dim;
    c[0] =  0;
    c[1] = -1;
    c[2] =  1;

    c = ini_c + 18 * dim;
    c[0] =  0;
    c[1] =  1;
    c[2] = -1;
    
    *one_over_c_s2 = 3.0;	// Inverso da vel. do som ao quadrado

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    W[0] = 1. / 3.;

    for ( int i = 1; i < 7; i++ ) W[i] = 1. / 18.;

    for ( int i = 7; i < nvel; i++ ) W[i] = 1. / 36.;       
}

//================================================================================================//




//===================== Define os vetores da rede D3Q15 ==========================================//
//
//      Input: pointer to the vectors
//      Output:
//
//================================================================================================//

void def_lattice_d3q15 ( double  *ini_c, double *W, double *one_over_c_s2 )
{
    double  *c;

    //------------- |c| = 0 ----------------------//

    c = ini_c;
    c[0] = 0;
    c[1] =  0;
    c[2] = 0;

    //------------- |c| = 1 ----------------------//

    c = ini_c + 1 * dim;
    c[0] =  1;
    c[1] =  0;
    c[2] =  0;

    c = ini_c + 2 * dim;
    c[0] = -1;
    c[1] =  0;
    c[2] =  0;

    c = ini_c + 3 * dim;
    c[0] =  0;
    c[1] =  1;
    c[2] =  0;

    c = ini_c + 4 * dim;
    c[0] =  0;
    c[1] = -1;
    c[2] =  0;

    c = ini_c + 5 * dim;
    c[0] =  0;
    c[1] =  0;
    c[2] =  1;

    c = ini_c + 6 * dim;
    c[0] =  0;
    c[1] =  0;
    c[2] = -1;

    //------------ |c| = sqrt(3) -----------------//

    c = ini_c + 7 * dim;
    c[0] =  1;
    c[1] =  1;
    c[2] =  1;

    c = ini_c + 8 * dim;
    c[0] = -1;
    c[1] = -1;
    c[2] = -1;

    c = ini_c + 9 * dim;
    c[0] = -1;
    c[1] =  1;
    c[2] =  1;

    c = ini_c + 10 * dim;
    c[0] =  1;
    c[1] = -1;
    c[2] = -1;

    c = ini_c + 11 * dim;
    c[0] =  1;
    c[1] = -1;
    c[2] =  1;

    c = ini_c + 12 * dim;
    c[0] = -1;
    c[1] =  1;
    c[2] = -1;

    c = ini_c + 13 * dim;
    c[0] =  1;
    c[1] =  1;
    c[2] = -1;

    c = ini_c + 14 * dim;
    c[0] = -1;
    c[1] = -1;
    c[2] =  1;
    
    *one_over_c_s2 = 3.0;	// Inverso da vel. do som ao quadrado 

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    W[0] = 2. / 9.;

    for ( int i = 1; i < 7; i++ ) W[i] = 1. / 9.;

    for ( int i = 7; i < nvel; i++ ) W[i] = 1. / 72.;

}

//================================================================================================//




//===================== Define os vetores da rede D3V59 ==========================================//
//
//      Input: pointer to the vectors
//      Output:
//
//================================================================================================//

void def_lattice_d3v59 ( double  *ini_c, double *W, double *one_over_c_s2 )
{
    double  *c;

    //------------- |c| = 0 ----------------------//

    c = ini_c;
    c[0] = 0;
    c[1] = 0;
    c[2] = 0;

    //------------- |c| = 1 ----------------------//

    c = ini_c + 1 * dim;
    c[0] = 1;
    c[1] = 0;
    c[2] = 0;

    c = ini_c + 2 * dim;
    c[0] = 0;
    c[1] = 1;
    c[2] = 0;

    c = ini_c + 3 * dim;
    c[0] = 0;
    c[1] = 0;
    c[2] = 1;

    c = ini_c + 4 * dim;
    c[0] = -1;
    c[1] = 0;
    c[2] = 0;

    c = ini_c + 5 * dim;
    c[0] = 0;
    c[1] = -1;
    c[2] = 0;

    c = ini_c + 6 * dim;
    c[0] = 0;
    c[1] = 0;
    c[2] = -1;

    //------------ |c| -> sqrt(2) -----------------//

    c = ini_c + 7 * dim;
    c[0] = 1;
    c[1] = 1;
    c[2] = 0;

    c = ini_c + 8 * dim;
    c[0] = -1;
    c[1] = 1;
    c[2] = 0;

    c = ini_c + 9 * dim;
    c[0] = -1;
    c[1] = -1;
    c[2] = 0;

    c = ini_c + 10 * dim;
    c[0] = 1;
    c[1] = -1;
    c[2] = 0;

    c = ini_c + 11 * dim;
    c[0] = 0;
    c[1] = 1;
    c[2] = -1;

    c = ini_c + 12 * dim;
    c[0] = 0;
    c[1] = 1;
    c[2] = 1;

    c = ini_c + 13 * dim;
    c[0] = 0;
    c[1] = -1;
    c[2] = 1;

    c = ini_c + 14 * dim;
    c[0] = 0;
    c[1] = -1;
    c[2] = -1;

    c = ini_c + 15 * dim;
    c[0] = 1;
    c[1] = 0;
    c[2] = -1;

    c = ini_c + 16 * dim;
    c[0] = -1;
    c[1] = 0;
    c[2] = -1;

    c = ini_c + 17 * dim;
    c[0] = -1;
    c[1] = 0;
    c[2] = 1;

    c = ini_c + 18 * dim;
    c[0] = 1;
    c[1] = 0;
    c[2] = 1;

    //------------ |c| -> sqrt(3) -----------------//

    c = ini_c + 19 * dim;
    c[0] = 1;
    c[1] = 1;
    c[2] = 1;

    c = ini_c + 20 * dim;
    c[0] = -1;
    c[1] = 1;
    c[2] = 1;

    c = ini_c + 21 * dim;
    c[0] = 1;
    c[1] = -1;
    c[2] = 1;

    c = ini_c + 22 * dim;
    c[0] = 1;
    c[1] = 1;
    c[2] = -1;

    c = ini_c + 23 * dim;
    c[0] = -1;
    c[1] = -1;
    c[2] = -1;

    c = ini_c + 24 * dim;
    c[0] = -1;
    c[1] = -1;
    c[2] = 1;

    c = ini_c + 25 * dim;
    c[0] = 1;
    c[1] = -1;
    c[2] = -1;

    c = ini_c + 26 * dim;
    c[0] = -1;
    c[1] = 1;
    c[2] = -1;

    //------------ |c| -> 2 -----------------//

    c = ini_c + 27 * dim;
    c[0] = 2;
    c[1] = 0;
    c[2] = 0;

    c = ini_c + 28 * dim;
    c[0] = 0;
    c[1] = 2;
    c[2] = 0;

    c = ini_c + 29 * dim;
    c[0] = 0;
    c[1] = 0;
    c[2] = 2;

    c = ini_c + 30 * dim;
    c[0] = -2;
    c[1] = 0;
    c[2] = 0;

    c = ini_c + 31 * dim;
    c[0] = 0;
    c[1] = -2;
    c[2] = 0;

    c = ini_c + 32 * dim;
    c[0] = 0;
    c[1] = 0;
    c[2] = -2;

    //------------ |c| -> 2 sqrt(2) -----------------//

    c = ini_c + 33 * dim;
    c[0] = 2;
    c[1] = 2;
    c[2] = 0;

    c = ini_c + 34 * dim;
    c[0] = -2;
    c[1] = 2;
    c[2] = 0;

    c = ini_c + 35 * dim;
    c[0] = -2;
    c[1] = -2;
    c[2] = 0;

    c = ini_c + 36 * dim;
    c[0] = 2;
    c[1] = -2;
    c[2] = 0;

    c = ini_c + 37 * dim;
    c[0] = 0;
    c[1] = 2;
    c[2] = -2;

    c = ini_c + 38 * dim;
    c[0] = 0;
    c[1] = 2;
    c[2] = 2;

    c = ini_c + 39 * dim;
    c[0] = 0;
    c[1] = -2;
    c[2] = 2;

    c = ini_c + 40 * dim;
    c[0] = 0;
    c[1] = -2;
    c[2] = -2;

    c = ini_c + 41 * dim;
    c[0] = 2;
    c[1] = 0;
    c[2] = -2;

    c = ini_c + 42 * dim;
    c[0] = -2;
    c[1] = 0;
    c[2] = -2;

    c = ini_c + 43 * dim;
    c[0] = -2;
    c[1] = 0;
    c[2] = 2;

    c = ini_c + 44 * dim;
    c[0] = 2;
    c[1] = 0;
    c[2] = 2;

    //------------ |c| -> 2 sqrt(3) -----------------//

    c = ini_c + 45 * dim;
    c[0] = 2;
    c[1] = 2;
    c[2] = 2;

    c = ini_c + 46 * dim;
    c[0] = -2;
    c[1] = 2;
    c[2] = 2;

    c = ini_c + 47 * dim;
    c[0] = 2;
    c[1] = -2;
    c[2] = 2;

    c = ini_c + 48 * dim;
    c[0] = 2;
    c[1] = 2;
    c[2] = -2;

    c = ini_c + 49 * dim;
    c[0] = -2;
    c[1] = -2;
    c[2] = -2;

    c = ini_c + 50 * dim;
    c[0] = -2;
    c[1] = -2;
    c[2] = 2;

    c = ini_c + 51 * dim;
    c[0] = 2;
    c[1] = -2;
    c[2] = -2;

    c = ini_c + 52 * dim;
    c[0] = -2;
    c[1] = 2;
    c[2] = -2;

    //------------ |c| -> 3--------------------------//

    c = ini_c + 53 * dim;
    c[0] = 3;
    c[1] = 0;
    c[2] = 0;

    c = ini_c + 54 * dim;
    c[0] = 0;
    c[1] = 3;
    c[2] = 0;

    c = ini_c + 55 * dim;
    c[0] = 0;
    c[1] = 0;
    c[2] = 3;

    c = ini_c + 56 * dim;
    c[0] = -3;
    c[1] = 0;
    c[2] = 0;

    c = ini_c + 57 * dim;
    c[0] = 0;
    c[1] = -3;
    c[2] = 0;

    c = ini_c + 58 * dim;
    c[0] = 0;
    c[1] = 0;
    c[2] = -3;

	*one_over_c_s2 = 1.0 / ( 0.831334581 * 0.831334581 ); // Inverso da vel. do som ao quadrado

    //-------------- Inicializa os pesos de acordo com a rede ------------------------------------//

    W[0] = 0.0958789162377528327290944502318849561751026273923599561095430491230444536551731987;

    for ( int i = 1; i < 7; i++ )
    {
        W[i] = 0.0731047082129148391094174556584327643133045030377292215552109509429216828659792144;
    }

    for ( int i = 7; i < 19; i++ )
    {
        W[i] = 0.0034658897109338004402496848234008098849685648112917249734187912080101109030889258;
    }

    for ( int i = 19; i < 27; i++ )
    {
        W[i] = 0.0366108082044515378737437706207046704098270120824206614444657745330939788438163106;
    }

    for ( int i = 27; i < 33; i++ )
    {
        W[i] = 0.0159235232232059553213542734282058102046784142530580659486045998739107023207608027;
    }

    for ( int i = 33; i < 45; i++ )
    {
        W[i] = 0.0025248084510509439390875492857286780439383835632091129974817362098589311760311626;
    }

    for ( int i = 45; i < 53; i++ )
    {
        W[i] = 0.0000726968662515158634643662368250684494442072833085637149780734785414288702883519;
    }

    for ( int i = 53; i < 59; i++ )
    {
        W[i] = 0.0007658794393468397060938785130819717826577889071787436568677554780749114270180559;
    }
}

//================================================================================================//




//===================== Calcula a distribuição de equilíbrio =====================================//
//
//      Input: velocities, density, lattice vectors, distribution function
//      Output: equilibrium distribution
//
//================================================================================================//

void dist_eq ( double vx, double vy, double vz, double rho,  double *ini_c, double *feq,
               double *W, double one_over_c_s2 )
{
    double *c;

    double v[dim];

    v[0] = vx;
    v[1] = vy;
    v[2] = vz;

    double vquad = ( vx * vx + vy * vy + vz * vz );

    double cv;

    feq[0] = W[0] * rho * ( 1.0 - 0.5 * vquad * one_over_c_s2 );

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        cv = prod_int ( c, v );

        feq[i] = W[i] * rho * ( 1.0 + cv * one_over_c_s2
                                + 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2
                                - 0.5 * vquad * one_over_c_s2 );
    }
}

//================================================================================================//




//===================== Calcula a distribuição de equilíbrio com termo de Exclusão ===============//
//
//      Input: velocities, density, lattice vectors, distribution function
//      Output: equilibrium distribution
//
//================================================================================================//

void dist_eq_E ( double vx, double vy, double vz, double rho, double gp_x, double gp_y, double gp_z,  
				double tau, double a, double *ini_c, double *feq, double *W, double one_over_c_s2 )
{
    double *c;

    double v[dim];

    v[0] = vx;
    v[1] = vy;
    v[2] = vz;
    
    double vquad = ( vx * vx + vy * vy + vz * vz );

    double cv, gpc, gpv;
    
    double gp[dim];
    
    gp[0] = gp_x;
    gp[1] = gp_y;
    gp[2] = gp_z;
    
    gpv = prod_int ( gp, v );

    feq[0] = W[0] * rho * ( 1.0 - 0.5 * vquad * one_over_c_s2 + tau * ( a * gpv ) );

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        cv = prod_int ( c, v );
        
        gpc = prod_int ( gp, c );
        
        feq[i] = W[i] * rho * ( 1.0 + cv * one_over_c_s2
                                + 0.5 * cv * cv * one_over_c_s2 * one_over_c_s2
                                - 0.5 * vquad * one_over_c_s2 
                                - tau * ( a*a*a * gpc * cv - a * gpv ) );
    }
}

//================================================================================================//




//===================== Calcula a distribuição de equilíbrio térmica (D3V59) =====================//
//
//      Input: velocities, density, theta, lattice vectors, distribution function
//      Output: equilibrium distribution
//
//================================================================================================//

void dist_eq_thermal ( double vx, double vy, double vz, double rho, double theta, double *ini_c,
                       double *feq, double *W )
{

    //====== ISOTÉRMICO =====//

    theta = 1.0 / 1.444;

    //=======================//

    double *c;

    double u[dim];

    u[0] = vx;
    u[1] = vy;
    u[2] = vz;

    double a = 1.2028851233102617112398875023476154612257535459964765342960063948674139162796177820;

    double a2 = a * a;
    double a4 = a2 * a2;
    double a6 = a2 * a2 * a2;
    double a8 = a4 * a4;
    double uu = prod_int( u, u );

    feq[0] = W[0] * rho * ( 1. - (1./2.)*(a2 * uu + 3 * (a2 * theta - 1)) +
                            (15./8.) * (a2 * theta - 1) * (a2 * theta - 1)
                            + (5./4.)*(a2 * theta - 1) * a2 * uu +
                            (1./8.) * a4 * uu * uu );

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        double cc = prod_int( c, c );
        double uc = prod_int( c, u );
        double cu = uc;

        feq[i] = W[i] * rho * ( 1. + a2 * cu + (1./2.)*( a4 * cu * cu  - a2 * uu
                                + (a2 * theta - 1) * (a2 * cc - 3))
                                + (1./6.) * a2 * cu * ( a4 * cu * cu - 3. * a2 * uu
                                        + 3. * (a2 * theta - 1) * (a2 * cc - 5. ))
                                + (15./8.) * (a2 * theta - 1) * (a2 * theta - 1)
                                + (5./4.)*(a2 * theta - 1) * a2 * uu
                                + (1./8.) * a4 * uu * uu
                                - (7./4.)*(a2 * theta - 1) * a4 * uc * uc
                                - (1./4.) * a6 * uu * uc * uc
                                - (5./4.) * (a2 * theta - 1)*(a2 * theta - 1) * a2 * cc
                                - (1./4.) * (a2 * theta
                                             - 1) * a4 * uu * cc + (1./4.)*(a2 * theta
                                                     - 1) * a6 * uc * uc * cc
                                + (1./28.) * a8 * uu * uc * uc * cc
                                + (1./8.)*(a2 * theta - 1)*(a2 * theta - 1) * a4 * cc * cc
                                - (1./280.) * a8 * uu * uu * cc * cc );
    }
}

//================================================================================================//





//===================== Calcula densidade e velocidades ==========================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density
//
//================================================================================================//

void calcula ( double *f, double *ini_c, double *vx, double *vy, double *vz, double *rho )
{
    double *c;

    double mx = 0.0;
    double my = 0.0;
    double mz = 0.0;

    double one_over_rho;
    
    double rho_tmp = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        rho_tmp = rho_tmp + f[i];

        mx = mx + c[0] * f[i];
        my = my + c[1] * f[i];
        mz = mz + c[2] * f[i];
    }

    if ( rho_tmp )
    {
        one_over_rho = 1.0 / rho_tmp;

        *vx = mx * one_over_rho;
        *vy = my * one_over_rho;
        *vz = mz * one_over_rho;
    }
    else
    {
        *vx = 0;
        *vy = 0;
        *vz = 0;
    }
    
    *rho = rho_tmp;
}

//================================================================================================//




//===================== Calcula o tensor c_alpha c_beta ==========================================//
//
//      Input: lattice vectors
//      Output: tensor cc
//
//================================================================================================//

void calcula_Q ( double* ini_c, double *ini_Q )
{
	double c_s2 = 1.0 / one_over_c_s2;

    for ( int i = 0 ; i < nvel; i++ )
    {
        double* c = ini_c + i * dim;

        for ( int alpha = 0; alpha < dim; alpha++ )
        {
            for ( int beta = 0; beta < dim; beta++ )
            {

                double* Q = ini_Q + alpha + beta * dim + i * dim * dim;

                *Q = c[alpha] * c[beta] - ( c_s2 ) * ( alpha == beta ) ;
            }
        }
    }
}

//================================================================================================//




//===================== Calcula densidade e velocidades ==========================================//
//
//      Input: distribution function, lattice vectors
//      Output: velocities, density
//
//================================================================================================//

void calcula_thermal ( double *f, double *ini_c, double *vx, double *vy, double *vz,
                       double *rho, double *theta )
{
    double *c;
    
    double vx_tmp = 0.0;
    double vy_tmp = 0.0;
    double vz_tmp = 0.0;

    double mx = 0.0;
    double my = 0.0;
    double mz = 0.0;

    double energia = 0.0;

    double one_over_rho;
    
    double rho_tmp = 0.0;

    for ( int i = 0 ; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        rho_tmp = rho_tmp + f[i];

        mx = mx + c[0] * f[i];
        my = my + c[1] * f[i];
        mz = mz + c[2] * f[i];
    }

    if ( rho_tmp )
    {
        one_over_rho = 1.0 / ( rho_tmp );

        vx_tmp = mx * one_over_rho;
        vy_tmp = my * one_over_rho;
        vz_tmp = mz * one_over_rho;

        for ( int i = 0 ; i < nvel; i++ )
        {
            c = ini_c + i * dim;

            double cx = c[0] - vx_tmp;
            double cy = c[1] - vy_tmp;
            double cz = c[2] - vz_tmp;

            energia = energia + f[i] * ( cx * cx + cy * cy + cz * cz );
        }

        energia = 0.5 * one_over_rho * energia;
    }
    else
    {
        vx_tmp = 0;
        vy_tmp = 0;
        vz_tmp = 0;
    }

    *theta = ( 2. / 3. ) * energia;
    
    *rho = rho_tmp;
    
    *vx = vx_tmp;
    *vy = vy_tmp;
    *vz = vz_tmp;
}

//================================================================================================//




//===================== Calcula a densidade no sítio =============================================//
//
//      Input: distribution function
//      Output: density
//
//================================================================================================//

double mass ( double *f )
{

    double rho = 0.0;

    for ( int i = 0 ; i < nvel; i++ ) rho = rho + f[i];

    return rho;
}

//================================================================================================//




//===================== Calcula o momento do sítio ===============================================//
//
//      Input: distribution function, lattice vectors
//      Output: momentum, *mx, *my, *mz
//
//================================================================================================//

void momentum ( double *f, double *ini_c, double *mx, double *my, double *mz )
{

    double *c;

    *mx = 0.0;
    *my = 0.0;
    *mz = 0.0;

    for ( int i = 1 ; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        *mx = *mx + c[0] * f[i];
        *my = *my + c[1] * f[i];
        *mz = *mz + c[2] * f[i];
    }
}

//================================================================================================//




//===================== Calcula um vetor unitário ================================================//
//
//      Input: vector
//      Output: unit vector, *ux, *uy, *uz
//
//================================================================================================//

void vet_unit ( double vx, double vy, double vz, double *ux, double *uy, double *uz )
{
    double modulo = ( sqrt ( vx * vx + vy * vy + vz * vz ) );

    double one_over_modulo;

    if ( modulo )
    {
        one_over_modulo = 1.0 / modulo;

        *ux = vx * one_over_modulo;
        *uy = vy * one_over_modulo;
        *uz = vz * one_over_modulo;
    }
    else
    {
        *ux = *uy = *uz = 0.0;
    }
}

//================================================================================================//




//===================== Calcula o operador BGK ===================================================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//================================================================================================//

void bgk_op ( double *f, double *ini_c, double vx, double vy, double vz, double rho,
              double acc_x, double acc_y, double acc_z, double tau, double *op_col,
              double *W, double one_over_c_s2 )
{
    double one_over_tau = 1.0 / tau;

    double *f_eq = new double[nvel];

    double vx_acc = vx + acc_x * tau;

    double vy_acc = vy + acc_y * tau;

    double vz_acc = vz + acc_z * tau;

    dist_eq ( vx_acc, vy_acc, vz_acc, rho, ini_c, f_eq, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        op_col[i] = ( f_eq[i] - f[i] ) * one_over_tau;
    }

	delete[] f_eq;
}

//================================================================================================//




//===================== Calcula o operador BGK simétrico (modelo TRT) ============================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//================================================================================================//

void bgk_even ( double *f, double *ini_c, double vx, double vy, double vz, double rho,
                double acc_x, double acc_y, double acc_z, double tau, double *op_col,
                double *W, double one_over_c_s2 )
{
    double f_eq[nvel];

    double vx_acc = vx + acc_x * tau;

    double vy_acc = vy + acc_y * tau;

    double vz_acc = vz + acc_z * tau;

    dist_eq ( vx_acc, vy_acc, vz_acc, rho, ini_c, f_eq, W, one_over_c_s2 );

    op_col[0] = ( f_eq[0] - f[0] ) / ( tau );

    for ( int i = 1; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] + f[i+1] ) / 2.0;

        double f_eq_mais = ( f_eq[i] + f_eq[i+1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) / tau;
    }

    for ( int i = 2; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] + f[i-1] ) / 2.0;

        double f_eq_mais = ( f_eq[i] + f_eq[i-1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) / tau;
    }
}

//================================================================================================//




//===================== Calcula o operador BGK simétrico (modelo TRT) - Exclusão por volume ======//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//================================================================================================//

void bgk_even_E ( double *f, double *ini_c, double vx, double vy, double vz, double rho,
					double acc_x, double acc_y, double acc_z, double gp_x, double gp_y, double gp_z, 
					double tau, double *op_col,	double *W, double one_over_c_s2 )
{
    double f_eq[nvel];

    double vx_acc = vx + acc_x * tau;

    double vy_acc = vy + acc_y * tau;

    double vz_acc = vz + acc_z * tau;
    
    double a = sqrt ( one_over_c_s2 );

    dist_eq_E ( vx_acc, vy_acc, vz_acc, rho, gp_x, gp_y, gp_z, tau, a, ini_c, f_eq, W, 
				one_over_c_s2 );

    op_col[0] = ( f_eq[0] - f[0] ) / ( tau );

    for ( int i = 1; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] + f[i+1] ) / 2.0;

        double f_eq_mais = ( f_eq[i] + f_eq[i+1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) / tau;
    }

    for ( int i = 2; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] + f[i-1] ) / 2.0;

        double f_eq_mais = ( f_eq[i] + f_eq[i-1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) / tau;
    }
}

//================================================================================================//




//===================== Calcula o operador BGK antisimétrico (modelo TRT) ========================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//================================================================================================//

void bgk_odd ( double *f, double *ini_c, double vx, double vy,
               double vz, double rho, double acc_x, double acc_y, double acc_z,
               double tau, double *op_col, double *W, double one_over_c_s2 )
{
    double f_eq[nvel];

    double vx_acc = vx + acc_x * tau;

    double vy_acc = vy + acc_y * tau;

    double vz_acc = vz + acc_z * tau;

    dist_eq ( vx_acc, vy_acc, vz_acc, rho, ini_c, f_eq, W, one_over_c_s2 );

    op_col[0] =  0.0;

    for ( int i = 1; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] - f[i+1] ) / 2.0;

        double f_eq_mais = ( f_eq[i] - f_eq[i+1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) / tau;
    }

    for ( int i = 2; i < nvel; i = i + 2 )
    {
        double f_mais =  ( f[i] - f[i-1] ) / 2.0;

        double f_eq_mais = ( f_eq[i] - f_eq[i-1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) / tau;
    }

}

//================================================================================================//




//===================== Calcula o operador BGK antisimétrico (modelo TRT) ========================//
//
//      Input: distribution function, lattice vectors, velocities, density, aceleration,
//              relaxation time
//      Output: BGK operator
//
//================================================================================================//

void bgk_odd_E ( double *f, double *ini_c, double vx, double vy, double vz, double rho, 
				double acc_x, double acc_y, double acc_z, double gp_x, double gp_y, 
				double gp_z, double tau, double *op_col, double *W, double one_over_c_s2 )
{
    double f_eq[nvel];

    double vx_acc = vx + acc_x * tau;

    double vy_acc = vy + acc_y * tau;

    double vz_acc = vz + acc_z * tau;
    
    double a = sqrt ( one_over_c_s2 );

    dist_eq_E ( vx_acc, vy_acc, vz_acc, rho, gp_x, gp_y, gp_z, tau, a, ini_c, f_eq, W, 
				one_over_c_s2 );

    op_col[0] =  0.0;

    for ( int i = 1; i < nvel; i = i + 2 )
    {
        double f_mais = ( f[i] - f[i+1] ) / 2.0;

        double f_eq_mais = ( f_eq[i] - f_eq[i+1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) / tau;
    }

    for ( int i = 2; i < nvel; i = i + 2 )
    {
        double f_mais =  ( f[i] - f[i-1] ) / 2.0;

        double f_eq_mais = ( f_eq[i] - f_eq[i-1] ) / 2.0;

        op_col[i] = ( f_eq_mais - f_mais ) / tau;
    }

}

//================================================================================================//




//===================== Etapa de colisão usando BGK ==============================================//
//
//      Input: distribution function, lattice vectors, acceleration in the x,y,z directions
//      Output: pos-collisional distribution function
//
//================================================================================================//

void bgk_collision ( double *f, double *ini_c, double tau, double acc_x, double acc_y,
                     double acc_z, double *W, double one_over_c_s2 )
{
    double *f_eq = new double[nvel];

    double *op_col = new double[nvel];

    double vx, vy, vz, ro;

    double one_over_tau = 1.0 / tau;

    calcula ( f, ini_c, &vx, &vy, &vz, &ro );

    double vx_alt = vx + acc_x * tau;

    double vy_alt = vy + acc_y * tau;

    double vz_alt = vz + acc_z * tau;

    dist_eq ( vx_alt, vy_alt, vz_alt, ro, ini_c, f_eq, W, one_over_c_s2 );

    for (int i=0; i<nvel; i++)
    {
        op_col[i] = ( f_eq[i] - f[i] ) * one_over_tau;

        f[i] = f[i] + op_col[i];
    }
    
	delete[] f_eq;
	
	delete[] op_col;
}

//================================================================================================//




//===================== Etapa de colisão com regularização - Latt & Chopard ======================//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//================================================================================================//


void regularized_collision ( double f[nvel], double *ini_c, double *ini_Q, double tau, double acc_x,
								double acc_y, double acc_z )
{
    const double c_s2 = 1. / one_over_c_s2;

    double omega = 1.0 / tau;

    double ux, uy, uz, rho;

    calcula ( f, ini_c, &ux, &uy, &uz, &rho );
    
    ux = ux + acc_x;	// Inclui aceleração
    
    ux = ux + acc_x;

    ux = ux + acc_x;    

    //----------- Calcula a distribuição de equilíbrio -------------------------------------------//

    double f_eq[nvel];

    dist_eq( ux, uy, uz, rho, ini_c, f_eq, W, one_over_c_s2);

    //----------- Calcula a distribuição de não-equilíbrio ---------------------------------------//

    double f_neq[nvel];

    for ( int i = 0; i < nvel; i++ ) f_neq[i] = f[i] - f_eq[i];

    //----------- Calcula o fluxo de não-equilíbio -----------------------------------------------//

    double Pi_xx = 0.0;

    double Pi_yy = 0.0;
    
    double Pi_zz = 0.0;

    double Pi_xy = 0.0;
    
    double Pi_xz = 0.0;
    
    double Pi_yz = 0.0;

    double *c = ini_c;

    for ( int i = 0; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        Pi_xx = Pi_xx + f_neq[i] * c[0] * c[0] - c_s2 * f_neq[i];

        Pi_yy = Pi_yy + f_neq[i] * c[1] * c[1] - c_s2 * f_neq[i];

		Pi_zz = Pi_zz + f_neq[i] * c[2] * c[2] - c_s2 * f_neq[i];

        Pi_xy = Pi_xy + f_neq[i] * c[0] * c[1];
        
        Pi_xz = Pi_xz + f_neq[i] * c[0] * c[2];
        
        Pi_yz = Pi_yz + f_neq[i] * c[1] * c[2];
    }

    //----------- Calcula f(1) -------------------------------------------------------------------//

    const double one_over_cs4 = one_over_c_s2 * one_over_c_s2;

    double f_1[nvel];

    for ( int i = 0; i < nvel; i ++ )
    {
        double* Q = ini_Q + i * dim * dim;

        f_1[i] = 0.5 * W[i] * one_over_cs4 * ( Q[0] * Pi_xx + Q[4] * Pi_yy + Q[8] * Pi_zz 
									+ 2 * Q[1] * Pi_xy + 2 * Q[2] * Pi_xz + 2 * Q[5] * Pi_yz );         
    }

    //-------------- Colisão ---------------------------------------------------------------------//

    for ( int i = 0; i < nvel; i++ )	f[i] = f_eq[i] + ( 1.0 - omega ) * f_1[i];

}

//================================================================================================//




//===================== Etapa de colisão com regularização - Latt & Chopard + TRT ================//
//
//      Input: distribution function, lattice vectors
//      Output: pos-collisional distribution function
//
//================================================================================================//


void regularized_TRT_collision ( double f[nvel], double *ini_c, double *ini_Q, double tau_even, 
									double tau_odd )
{
    const double c_s2 = 1. / one_over_c_s2;

    double omega_even = 1.0 / tau_even;
    
    double omega_odd = 1.0 / tau_odd;

    double ux, uy, uz, rho;

    calcula ( f, ini_c, &ux, &uy, &uz, &rho );

    //----------- Calcula a distribuição de equilíbrio -------------------------------------------//

    double f_eq[nvel];

    dist_eq( ux, uy, uz, rho, ini_c, f_eq, W, one_over_c_s2);
    
    //----------- Calcula as distribuições simétrica e anti-simétrica ----------------------------//
    
    double f_plus[nvel];
    
    double f_minus[nvel];
    
    double f_eq_plus[nvel];
    
    double f_eq_minus[nvel];
    
    f_plus[0] = f[0];
    
    f_minus[0] = 0.0;
    
    f_eq_plus[0] = f_eq[0];
    
    f_eq_minus[0] = 0.0;
    
    for ( int i = 1; i < nvel; i = i + 2 )
    {
		f_plus[i] = 0.5 * ( f[i] + f[i+1] );
		
		f_minus[i] = 0.5 * ( f[i] - f[i+1] );
		
		f_eq_plus[i] = 0.5 * ( f_eq[i] + f_eq[i+1] );
		
		f_eq_minus[i] = 0.5 * ( f_eq[i] - f_eq[i+1] );
	} 
	
	for ( int i = 2; i < nvel; i = i + 2 )
    {
		f_plus[i] = 0.5 * ( f[i] + f[i-1] );
		
		f_minus[i] = 0.5 * ( f[i] - f[i-1] );
		
		f_eq_plus[i] = 0.5 * ( f_eq[i] + f_eq[i-1] );
		
		f_eq_minus[i] = 0.5 * ( f_eq[i] - f_eq[i-1] ); 
	}
    
    //----------- Calcula a distribuição de não-equilíbrio ---------------------------------------//

    double f_neq_plus[nvel];
    
    double f_neq_minus[nvel];

    for ( int i = 0; i < nvel; i++ )
    {
		f_neq_plus[i] = f_plus[i] - f_eq_plus[i];
		
		f_neq_minus[i] = f_minus[i] - f_eq_minus[i];
	}

    //----------- Calcula o fluxo de não-equilíbio simétrico -------------------------------------//

    double Pi_xx_plus = 0.0;
    double Pi_yy_plus = 0.0;    
    double Pi_zz_plus = 0.0;
    
    double Pi_xy_plus = 0.0;    
    double Pi_xz_plus = 0.0;    
    double Pi_yz_plus = 0.0;

    double *c = ini_c;

    for ( int i = 0; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        Pi_xx_plus = Pi_xx_plus + f_neq_plus[i] * c[0] * c[0] - c_s2 * f_neq_plus[i];
        Pi_yy_plus = Pi_yy_plus + f_neq_plus[i] * c[1] * c[1] - c_s2 * f_neq_plus[i];
		Pi_zz_plus = Pi_zz_plus + f_neq_plus[i] * c[2] * c[2] - c_s2 * f_neq_plus[i];

        Pi_xy_plus = Pi_xy_plus + f_neq_plus[i] * c[0] * c[1];        
        Pi_xz_plus = Pi_xz_plus + f_neq_plus[i] * c[0] * c[2];        
        Pi_yz_plus = Pi_yz_plus + f_neq_plus[i] * c[1] * c[2];
    }
    
    //----------- Calcula o fluxo de não-equilíbio anti-simétrico --------------------------------//

    double Pi_xx_minus = 0.0;
    double Pi_yy_minus = 0.0;    
    double Pi_zz_minus = 0.0;
    
    double Pi_xy_minus = 0.0;    
    double Pi_xz_minus = 0.0;    
    double Pi_yz_minus = 0.0;

    for ( int i = 0; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        Pi_xx_minus = Pi_xx_minus + f_neq_minus[i] * c[0] * c[0] - c_s2 * f_neq_minus[i];
        Pi_yy_minus = Pi_yy_minus + f_neq_minus[i] * c[1] * c[1] - c_s2 * f_neq_minus[i];
		Pi_zz_minus = Pi_zz_minus + f_neq_minus[i] * c[2] * c[2] - c_s2 * f_neq_minus[i];

        Pi_xy_minus = Pi_xy_minus + f_neq_minus[i] * c[0] * c[1];        
        Pi_xz_minus = Pi_xz_minus + f_neq_minus[i] * c[0] * c[2];        
        Pi_yz_minus = Pi_yz_minus + f_neq_minus[i] * c[1] * c[2];
    }
    
    //----------- Calcula f(1) simétrico ---------------------------------------------------------//

    const double one_over_cs4 = one_over_c_s2 * one_over_c_s2;

    double f_1_plus[nvel];

    for ( int i = 0; i < nvel; i ++ )
    {
        double* Q = ini_Q + i * dim * dim;

        f_1_plus[i] = 0.5 * W[i] * one_over_cs4 * ( Q[0] * Pi_xx_plus + Q[4] * Pi_yy_plus 
						+ Q[8] * Pi_zz_plus + 2 * Q[1] * Pi_xy_plus + 2 * Q[2] * Pi_xz_plus 
						+ 2 * Q[5] * Pi_yz_plus );         
    }
    
    //----------- Calcula f(1) anti-simétrico ----------------------------------------------------//

    double f_1_minus[nvel];

    for ( int i = 0; i < nvel; i ++ )
    {
        double* Q = ini_Q + i * dim * dim;

        f_1_minus[i] = 0.5 * W[i] * one_over_cs4 * ( Q[0] * Pi_xx_minus + Q[4] * Pi_yy_minus 
						+ Q[8] * Pi_zz_minus + 2 * Q[1] * Pi_xy_minus + 2 * Q[2] * Pi_xz_minus 
						+ 2 * Q[5] * Pi_yz_minus );         
    }

    //-------------- Colisão ---------------------------------------------------------------------//

    for ( int i = 0; i < nvel; i++ )	f[i] = f_eq[i] + ( 1.0 - omega_even ) * f_1_plus[i]
											    + ( 1.0 - omega_odd ) * f_1_minus[i];

}

//================================================================================================//




//===================== Etapa de colisão usando BGK para rede d3v59 ==============================//
//
//      Input: distribution function, lattice vectors, acceleration in the x,y,z directions
//      Output: pos-collisional distribution function
//
//================================================================================================//

void bgk_collision_thermal ( double *f, double *ini_c, double tau, double acc_x, double acc_y,
                             double acc_z, double *W, double one_over_c_s2 )
{
    double *f_eq = new double[nvel];

    double *op_col = new double[nvel];

    double vx, vy, vz, rho, theta;

    double one_over_tau = 1.0 / tau;

    calcula_thermal ( f, ini_c, &vx, &vy, &vz, &rho, &theta );

    double vx_alt = vx + acc_x * tau;

    double vy_alt = vy + acc_y * tau;

    double vz_alt = vz + acc_z * tau;

    dist_eq_thermal ( vx_alt, vy_alt, vz_alt, rho, theta, ini_c, f_eq, W );

    for (int i=0; i<nvel; i++)
    {
        op_col[i] = ( f_eq[i] - f[i] ) * one_over_tau;

        f[i] = f[i] + op_col[i];
    }

    delete[] f_eq;
    
    delete[] op_col;
}

//================================================================================================//




//===================== Etapa de colisão usando TRT ==============================================//
//
//      Input: distribution function, lattice vectors, relaxation times
//      Output: pos-collisional distribution function
//
//================================================================================================//

void trt_collision ( double *f, double *ini_c, double tau_sim, double tau_ant, double acc_x,
                     double acc_y, double acc_z, double *W, double one_over_c_s2 )
{
    double vx, vy, vz, ro;

    calcula ( f, ini_c, &vx, &vy, &vz, &ro );

    double *op_col_sim = new double[nvel];

    bgk_even ( f, ini_c, vx, vy, vz, ro, acc_x, acc_y, acc_z, tau_sim, op_col_sim, W,
               one_over_c_s2 );

    double *op_col_ant = new double[nvel];

    bgk_odd ( f, ini_c, vx, vy, vz, ro, acc_x, acc_y, acc_z, tau_ant, op_col_ant, W,
              one_over_c_s2 );

	
    for ( int i = 0; i < nvel; i++ )
    {
        f[i] = f[i] + op_col_sim[i] + op_col_ant[i];
    }
    delete[] op_col_sim;
    
    delete[] op_col_ant;
}

//================================================================================================//




//===================== Etapa de colisão usando TRT (para advecção/difusão) ======================//
//
//      Input: distribution function, lattice vectors, relaxation times, velocities
//      Output: pos-collisional distribution function
//
//================================================================================================//

void trt_collision_theta ( double *g, double *ini_c, double tau_sim, double tau_ant, double vx,
                     double vy, double vz, double *W, double one_over_c_s2 )
{
	double acc_x = 0.0;
	double acc_y = 0.0;
	double acc_z = 0.0;
	
    double theta = mass( g );

    double *op_col_sim = new double[nvel];

    bgk_even ( g, ini_c, vx, vy, vz, theta, acc_x, acc_y, acc_z, tau_sim, op_col_sim, W,
               one_over_c_s2 );

    double *op_col_ant = new double[nvel];

    bgk_odd ( g, ini_c, vx, vy, vz, theta, acc_x, acc_y, acc_z, tau_ant, op_col_ant, W,
              one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        g[i] = g[i] + op_col_sim[i] + op_col_ant[i];

    }
    
    delete[] op_col_sim;
    
    delete[] op_col_ant;
}

//================================================================================================//




//===================== Determina se há variação na posição da interface =========================//
//
//      Input: geometri, mass fractions, old mass fractions, nx, ny, nz
//      Output: variação da posição da interface
//
//================================================================================================//

double var_inter ( int *ini_meio, double *ini_omega_RB, double *ini_omega_RB_old, int nx, 
							int ny, int nz )
{
	double gt_delta = 0.0;
	
	for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int *meio = ini_meio + x + y * nx + z * nx * ny;

                if ( *meio )
                {
					double *omega_RB = ini_omega_RB + x  + y * nx + z * nx * ny ;
					
					double *omega_RB_old = ini_omega_RB_old + x  + y * nx + z * nx * ny ;
					
					double delta = fabs ( fabs( *omega_RB ) - fabs( *omega_RB_old ) );
					
					if ( delta > gt_delta ) gt_delta = delta;
				}
			}
		}
	}
	
	return gt_delta;
}

//================================================================================================//




//===================== Encontra as direções opostas aos vet. da rede ============================//
//
//      Input: lattice vectors
//      Output: 
//
//================================================================================================//

void i_opposit ( double *ini_c, int *i_op )
{
	for ( int i = 1; i < nvel; i++ )
    {
        double *c = ini_c + i * dim;

        int cx_i = sum_round ( c[0] );
        int cy_i = sum_round ( c[1] );
        int cz_i = sum_round ( c[2] );

        for ( int j = 1; j < nvel; j++ )
        {
            c = ini_c + j * dim;

            int cx_j = sum_round ( c[0] );
            int cy_j = sum_round ( c[1] );
            int cz_j = sum_round ( c[2] );

            if ( cx_j == -cx_i  && cy_j == -cy_i && cz_j == -cz_i )
            {
                i_op[i] = j;
            }
        }
    }
}

//================================================================================================//




//===================== Define os sítios para a etapa de propagação ==============================//
//
//      Input: lattice vectors, geometry, dimensions
//      Output: addresses, *ini_dir
//
//================================================================================================//

void def_dir_prop ( double *ini_c, int *ini_meio, const int *ini_dir, int nx, int ny, int nz )
{
    int i;
    double  *c;
    int *meio_local;
    int *meio_prop;

    int *dir;

    //---------------- Define os passos para propagação ------------------------------------------//
    
    double *step_x = new double[nvel];
    double *step_y = new double[nvel];
    double *step_z = new double[nvel];

    int *steps = new int[nvel];

    for ( i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        double  cx = c[0];
        double  cy = c[1];
        double  cz = c[2];

        int cx_int =  sum_round ( cx );
        int cy_int =  sum_round ( cy );
        int cz_int =  sum_round ( cz );

        steps[i] = abs ( cx_int );

        if ( abs ( cy_int ) > steps[i] )
        {
            steps[i] = abs ( cy_int );
        }

        if ( abs ( cz_int ) > steps[i] )
        {
            steps[i] = abs ( cz_int );
        }

        step_x[i] = cx / steps[i];
        step_y[i] = cy / steps[i];
        step_z[i] = cz / steps[i];
    }

    //-------------- Encontra as direções contrárias ---------------------------------------------//

    int *i_op = new int[nvel];
    
    i_opposit ( ini_c, i_op );

    //--------------------------------------------------------------------------------------------//

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio_local = ini_meio + x + y * nx + z * nx * ny;

                if ( *meio_local != 0 )
                {
                    dir = ( int* ) ini_dir + ( abs ( *meio_local ) - 1 ) * nvel;

                    meio_prop = meio_local;

                    dir[0] = ( abs( *meio_prop ) - 1 ) * nvel;

                    //----------------------------------------------------------------------------//

                    for ( i = 1; i < nvel; i++ )
                    {
                        int inv = 0; // número de inversões

                        double  stpx = step_x[i];
                        double  stpy = step_y[i];
                        double  stpz = step_z[i];

                        int x_prop = x;
                        int y_prop = y;
                        int z_prop = z;

                        double  x_prop_f = ( double ) x;
                        double  y_prop_f = ( double ) y;
                        double  z_prop_f = ( double ) z;

                        for ( int stp = 1; stp < steps[i] + 1; stp++ )
                        {
                            x_prop_f = x_prop_f + stpx;
                            y_prop_f = y_prop_f + stpy;
                            z_prop_f = z_prop_f + stpz;

                            x_prop = ( sum_round ( x_prop_f ) + nx ) % nx;
                            y_prop = ( sum_round ( y_prop_f ) + ny ) % ny;
                            z_prop = ( sum_round ( z_prop_f ) + nz ) % nz;

                            meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

                            if ( *meio_prop == 0 )
                            {
                                stpx = -stpx;
                                stpy = -stpy;
                                stpz = -stpz;

                                x_prop_f = x_prop_f + stpx;
                                y_prop_f = y_prop_f + stpy;
                                z_prop_f = z_prop_f + stpz;

                                x_prop = ( sum_round ( x_prop_f ) + nx ) % nx;
                                y_prop = ( sum_round ( y_prop_f ) + ny ) % ny;
                                z_prop = ( sum_round ( z_prop_f ) + nz ) % nz;

                                meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

                                inv++;
                            }
                        }

                        if ( inv % 2 == 0 )
                        {
                            dir[i] = i + ( abs ( *meio_prop ) - 1 ) * nvel;
                        }
                        else
                        {
                            dir[i] = i_op[i] + ( abs ( *meio_prop ) - 1 ) * nvel;
                        }
                    }

                    //----------------------------------------------------------------------------//
                }
            }
        }
    }
    
    delete[] step_x;
    delete[] step_y;
    delete[] step_z;
    
    delete[] steps;
    delete[] i_op;    
}

//================================================================================================//




//===================== Define os sítios para a etapa de propagação (com mediadores) =============//
//
//      Input: lattice vectors, geometry, dimensions
//      Output: addresses, *ini_dir, and momentum lost, *ini_qlost
//
//================================================================================================//

void def_dir_prop_med ( double *ini_c, int *ini_meio, const int *ini_dir, const bool* ini_solid, 
						int nx, int ny, int nz )
{
    int i;
    double  *c;
    int *meio_local;
    int *meio_prop;

    int *dir;
    bool* solid;

    //---------------- Define os passos para propagação ------------------------------------------//
    
    double *step_x = new double[nvel];
    double *step_y = new double[nvel];
    double *step_z = new double[nvel];

    int *steps = new int[nvel];

    for ( i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        double  cx = c[0];
        double  cy = c[1];
        double  cz = c[2];

        int cx_int =  sum_round ( cx );
        int cy_int =  sum_round ( cy );
        int cz_int =  sum_round ( cz );

        steps[i] = abs ( cx_int );

        if ( abs ( cy_int ) > steps[i] )
        {
            steps[i] = abs ( cy_int );
        }

        if ( abs ( cz_int ) > steps[i] )
        {
            steps[i] = abs ( cz_int );
        }

        step_x[i] = cx / steps[i];
        step_y[i] = cy / steps[i];
        step_z[i] = cz / steps[i];
    }

    //-------------- Encontra as direções contrárias ---------------------------------------------//

    int *i_op = new int[nvel];
    
    i_opposit ( ini_c, i_op );

    //--------------------------------------------------------------------------------------------//

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio_local = ini_meio + x + y * nx + z * nx * ny;

                if ( *meio_local )
                {
                    dir = ( int* ) ini_dir + ( *meio_local - 1 ) * nvel;
                    
                    solid = ( bool* ) ini_solid + ( *meio_local - 1 ) * nvel;

                    meio_prop = meio_local;

                    dir[0] = ( *meio_prop - 1 ) * nvel;
                    
                    solid[0] = 0;

                    //----------------------------------------------------------------------------//

                    for ( i = 1; i < nvel; i++ )
                    {
                        int inv = 0; // número de inversões

                        double  stpx = step_x[i];
                        double  stpy = step_y[i];
                        double  stpz = step_z[i];

                        int x_prop = x;
                        int y_prop = y;
                        int z_prop = z;

                        double  x_prop_f = ( double ) x;
                        double  y_prop_f = ( double ) y;
                        double  z_prop_f = ( double ) z;

                        for ( int stp = 1; stp < steps[i] + 1; stp++ )
                        {
                            x_prop_f = x_prop_f + stpx;
                            y_prop_f = y_prop_f + stpy;
                            z_prop_f = z_prop_f + stpz;

                            x_prop = ( sum_round ( x_prop_f ) + nx ) % nx;
                            y_prop = ( sum_round ( y_prop_f ) + ny ) % ny;
                            z_prop = ( sum_round ( z_prop_f ) + nz ) % nz;

                            meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

                            if ( *meio_prop == 0 )
                            {
                                stpx = -stpx;
                                stpy = -stpy;
                                stpz = -stpz;

                                x_prop_f = x_prop_f + stpx;
                                y_prop_f = y_prop_f + stpy;
                                z_prop_f = z_prop_f + stpz;

                                x_prop = ( sum_round ( x_prop_f ) + nx ) % nx;
                                y_prop = ( sum_round ( y_prop_f ) + ny ) % ny;
                                z_prop = ( sum_round ( z_prop_f ) + nz ) % nz;

                                meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

                                inv++;
                            }
                        }

                        if ( inv % 2 == 0 )
                        {
                            dir[i] = i + ( *meio_prop - 1 ) * nvel;
                        }
                        else
                        {
                            dir[i] = i_op[i] + ( *meio_prop - 1 ) * nvel;
                            
                            solid[i] = 1;
                        }
                    }

                    //----------------------------------------------------------------------------//
                }
            }
        }
    }
    
    delete[] step_x;
    delete[] step_y;
    delete[] step_z;
    
    delete[] steps;
    delete[] i_op;    
}

//================================================================================================//




//===================== Define os sítios para a etapa de propagação otimizada ====================//
//
//      Input: lattice vectors, geometry, dimensions
//      Output: addresses, *ini_dir
//
//================================================================================================//

void def_dir_prop_opt ( double *ini_c, int *ini_meio, int *ini_dir, int nx, int ny, int nz )
{
    int i;
    double  *c;
    int *meio_local;
    int *meio_prop;

    int *dir;

    //---------------- Define os passos para propagação ------------------------------------------//
    
    double *step_x = new double[nvel];
    double *step_y = new double[nvel];
    double *step_z = new double[nvel];

    int *steps = new int[nvel];

    for ( i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        double  cx = c[0];
        double  cy = c[1];
        double  cz = c[2];

        int cx_int =  sum_round ( cx );
        int cy_int =  sum_round ( cy );
        int cz_int =  sum_round ( cz );

        steps[i] = abs ( cx_int );

        if ( abs ( cy_int ) > steps[i] )
        {
            steps[i] = abs ( cy_int );
        }

        if ( abs ( cz_int ) > steps[i] )
        {
            steps[i] = abs ( cz_int );
        }

        step_x[i] = cx / steps[i];
        step_y[i] = cy / steps[i];
        step_z[i] = cz / steps[i];
    }

    //-------------- Encontra as direções contrárias ---------------------------------------------//

    int *i_op = new int[nvel];
    
    i_opposit ( ini_c, i_op );

    //--------------------------------------------------------------------------------------------//

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio_local = ini_meio + x + y * nx + z * nx * ny;

                if ( *meio_local )
                {
                    dir = ini_dir + ( *meio_local - 1 ) * nvel;

                    meio_prop = meio_local;

                    dir[0] = ( *meio_prop - 1 ) * nvel;

                    //----------------------------------------------------------------------------//

                    for ( i = 1; i < nvel; i++ )
                    {
                        int inv = 0; // número de inversões

                        double  stpx = step_x[i];
                        double  stpy = step_y[i];
                        double  stpz = step_z[i];

                        int x_prop = x;
                        int y_prop = y;
                        int z_prop = z;

                        double  x_prop_f = ( double ) x;
                        double  y_prop_f = ( double ) y;
                        double  z_prop_f = ( double ) z;

                        for ( int stp = 1; stp < steps[i] + 1; stp++ )
                        {
                            x_prop_f = x_prop_f + stpx;
                            y_prop_f = y_prop_f + stpy;
                            z_prop_f = z_prop_f + stpz;

                            x_prop = ( sum_round ( x_prop_f ) + nx ) % nx;
                            y_prop = ( sum_round ( y_prop_f ) + ny ) % ny;
                            z_prop = ( sum_round ( z_prop_f ) + nz ) % nz;

                            meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

                            if ( *meio_prop == 0 )
                            {
                                stpx = -stpx;
                                stpy = -stpy;
                                stpz = -stpz;

                                x_prop_f = x_prop_f + stpx;
                                y_prop_f = y_prop_f + stpy;
                                z_prop_f = z_prop_f + stpz;

                                x_prop = ( sum_round ( x_prop_f ) + nx ) % nx;
                                y_prop = ( sum_round ( y_prop_f ) + ny ) % ny;
                                z_prop = ( sum_round ( z_prop_f ) + nz ) % nz;

                                meio_prop = ini_meio + x_prop + y_prop * nx + z_prop * nx * ny;

                                inv++;
                            }
                        }

                        if ( inv % 2 == 0 )
                        {
                            dir[i_op[i]] = i + ( *meio_prop - 1 ) * nvel;
                        }
                        else
                        {
                            dir[i_op[i]] = i_op[i] + ( *meio_prop - 1 ) * nvel;
                        }
                    }

                    //----------------------------------------------------------------------------//
                }
            }
        }
    }
    
    delete[] step_x;
    delete[] step_y;
    delete[] step_z;
    
    delete[] steps;
    delete[] i_op;    
}

//================================================================================================//




//===================== Etapa de propagação para a distribuição de partículas ====================//
//
//      Input: distribution function (*ini_f), new distribution function (*ini_f_new),
//              geometry (*ini_meio), number of points
//      Output: updated distribution function (*ini_f_new)
//
//================================================================================================//

void propag_part ( int *ini_dir, double *ini_f, double *ini_f_new, int *ini_meio, int ptos_meio )
{
#pragma omp parallel for

    for ( int pto = 0; pto < ptos_meio; pto++ )
    {
        //--------------- Aponta os ponteiros ----------------------------------------------------//

        double *f = ini_f + ( pto ) * nvel;

        int *dir = ini_dir + ( pto ) * nvel;

        //---------------------------------------------------------------------------------------//

        for ( int i = 0; i < nvel; i++ )
        {
            double *f_new;

            f_new = ini_f_new + dir[i];

            *f_new = f[i];
        }
    }
}

//================================================================================================//




//===================== Etapa de propagação otimizada para a distribuição de partículas ==========//
//
//      Input: distribution function (*ini_f), new distribution function (*ini_f_new),
//              geometry (*ini_meio), number of points
//      Output: updated distribution function (*ini_f_new)
//
//================================================================================================//

void propag_part_opt ( int *ini_dir, double *ini_f, int *ini_meio, int ptos_meio )
{
#pragma omp parallel for

	for ( int i = 1; i < nvel; i = i + 2 )
    {
		for ( int pto = 0; pto < ptos_meio; pto++ )
		{
			int point = ( pto ) * nvel;
			
			//--------------- Aponta os ponteiros ------------------------------------------------//

			double *f_one = ini_f + point + i;

			int *dir = ini_dir + point;
			
			double *f_two = ini_f + dir[i];

			//------------------------------------------------------------------------------------//
        
            double f_swap = *f_one;

			*f_one = *f_two;
			
			*f_two = f_swap;
        }
    }
}

//================================================================================================//





//===================== Etapa de propagação para a um único sítio ================================//
//
//      Input: distribution function (*ini_f), new distribution function (*ini_f_new),
//              geometry (*ini_meio), point
//      Output: updated distribution function (*ini_f_new)
//
//================================================================================================//

void propag_part_site ( const int *ini_dir, double *ini_f, double *ini_f_new, int pto )
{

    //--------------- Aponta os ponteiros --------------------------------------------------------//

    double *f = ini_f + ( pto ) * nvel;

    int *dir = ( int* ) ini_dir + ( pto ) * nvel;

    //--------------------------------------------------------------------------------------------//
    
    for ( int i = 0; i < nvel; i++ )
    {
        double *f_new;

        f_new = ini_f_new + dir[i];

        *f_new = f[i];
    }
}

//================================================================================================//




//===================== Etapa de propagação para a um único sítio ================================//
//
//      Input: distribution function (*ini_f), new distribution function (*ini_f_new),
//              geometry (*ini_meio), point
//      Output: updated distribution function (*ini_f_new)
//
//================================================================================================//

void propag_g_site ( const int *ini_dir, double *ini_f, double *ini_f_new, const bool *ini_solid, 
					 double theta_wall, int pto, double* W )
{

    //--------------- Aponta os ponteiros --------------------------------------------------------//

    double *f = ini_f + ( pto ) * nvel;

    int *dir = ( int* ) ini_dir + ( pto ) * nvel;
    
    bool* solid = ( bool* ) ini_solid + ( pto ) * nvel;

    //--------------------------------------------------------------------------------------------//

    for ( int i = 0; i < nvel; i++ )
    {
        double *f_new;

        f_new = ini_f_new + dir[i];

        if ( solid[i] )	*f_new = -f[i] + 2.0 * W[i] * theta_wall;  

		else *f_new = f[i];
    }
}

//================================================================================================//




//===================== Etapa de propagação para a um único sítio (mediadores) ===================//
//
//      Input: distribution function (*ini_f), new distribution function (*ini_f_new),
//              geometry (*ini_meio), point
//      Output: updated distribution function (*ini_f_new)
//
//================================================================================================//

void propag_med_site ( int* ini_dir, double* ini_f, double* ini_f_new, bool* ini_solid, 
						double molhab, int pto, double* W )
{

    //--------------- Aponta os ponteiros --------------------------------------------------------//

    double* f = ini_f + ( pto ) * nvel;

    int* dir = ini_dir + ( pto ) * nvel;
    
    bool* solid = ini_solid + ( pto ) * nvel;

    //--------------------------------------------------------------------------------------------//

    for ( int i = 0; i < nvel; i++ )
    {
        double *f_new;

        f_new = ini_f_new + dir[i];
        
        if ( solid[i] )	*f_new = W[i] * molhab;  

		else *f_new = f[i];
    }
}

//================================================================================================//




//===================== Determina a posição da interface - y =====================================//
//
//      Input: geometry, distributions functions, lattice vectors, dimensions
//      Output: interface position
//
//================================================================================================//

double pos_interf_y ( int *ini_meio, double *ini_R,  double *ini_B, int x_0, int y_0, int y_f, 
						int z_0, int nx, int ny, int nz )
{
 
	double pos_y = 0.0;
	
	double sum_pond = 0.0;
	
	for ( int y = y_0; y < y_f; y++ )
    {
		int *meio = ini_meio + x_0 + y * nx + z_0 * ny * nx;
                
		if ( *meio )
		{
			double *R = ini_R + ( *meio - 1 ) * nvel;
			
			double *B = ini_B + ( *meio - 1 ) * nvel;

			double rho_R = mass( R );
            
            double rho_B = mass( B );
            
            double rho = rho_R + rho_B;

			double conc_R = rho_R / rho;
			
			double conc_B = 1.0 - conc_R;
			
			pos_y = pos_y + ( double ) y * conc_R * conc_B;
			
			sum_pond = sum_pond + conc_R * conc_B;
        }
    }  
    
    pos_y = pos_y / sum_pond;
    
    return pos_y;  
}

//================================================================================================//




//===================== Determina a posição da interface - x =====================================//
//
//      Input: geometry, distributions functions, lattice vectors, dimensions
//      Output: interface position
//
//================================================================================================//

double pos_interf_x ( int *ini_meio, double *ini_R,  double *ini_B, int x_0, int x_f, int y_0, 
						int z_0, int nx, int ny, int nz )
{
 
	double pos_x = 0.0;
	
	double sum_pond = 0.0;
	
	for ( int x = x_0; x < x_f; x++ )
    {
		int *meio = ini_meio + x + y_0 * nx + z_0 * ny * nx;
                
		if ( *meio )
		{
			double *R = ini_R + ( *meio - 1 ) * nvel;
			
			double *B = ini_B + ( *meio - 1 ) * nvel;

			double rho_R = mass( R );
            
            double rho_B = mass( B );
            
            double rho = rho_R + rho_B;

			double conc_R = rho_R / rho;
			
			double conc_B = 1.0 - conc_R;
			
			pos_x = pos_x + ( double ) x * conc_R * conc_B;
			
			sum_pond = sum_pond + conc_R * conc_B;
        }
    }  
    
    pos_x = pos_x / sum_pond;
    
    return pos_x;  
}

//================================================================================================//




//===================== Inverte a direção y (em um plano xz) =====================================//
//
//      Input: distribution function (*ini_f), geometry (*ini_meio), lattice vectors (*ini_c)
//              position of the plane xz, dimensions
//      Output:
//
//================================================================================================//

void inv_y ( double *ini_f, int *ini_meio, double *ini_c, int y, int nx, int ny, int nz )
{

    //-------------- Determina as direções com y invertido ---------------------------------------//

    int *i_yinv = new int[nvel];

    for ( int i = 1; i < nvel; i++ )
    {
        i_yinv[i] = 0;

        for ( int j = i; j < nvel; j++ )
        {
            double *c_i = ini_c + i * dim;

            double *c_j = ini_c + j * dim;

            if ( c_i[1] != 0 && c_i[1] == - c_j[1] && c_i[0] == c_j[0] && c_i[2] == c_j[2] )
            {
                i_yinv[i] = j;
            }
        }
    }

    //--------------------------------------------------------------------------------------------//

    for ( int z = 0; z < nz; z++ )
    {
        for ( int x = 0; x < nx; x++ )
        {
            int *meio = ini_meio + x + y * nx + z * nx * ny;

            if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                for ( int i = 1; i < nvel; i++ )
                {
                    if ( i_yinv[i] != 0 )
                    {
                        double cup = f[i];

                        f[i]  = f[ i_yinv[i] ];

                        f[ i_yinv[i] ] = cup;
                    }
                }
            }
        }
    }
    
    delete[] i_yinv;
}

//================================================================================================//




//===================== Inverte a direção z (em um plano xy) =====================================//
//
//      Input: distribution function (*ini_f), geometry (*ini_meio), lattice vectors (*ini_c)
//              position of the plane xz, dimensions
//      Output:
//
//================================================================================================//

void inv_z ( double *ini_f, int *ini_meio, double *ini_c, int z, int nx, int ny, int nz )
{

    //-------------- Determina as direções com z invertido ---------------------------------------//

    int *i_zinv = new int[nvel];

    for ( int i = 1; i < nvel; i++ )
    {
        i_zinv[i] = 0;

        for ( int j = i; j < nvel; j++ )
        {
            double *c_i = ini_c + i * dim;

            double *c_j = ini_c + j * dim;

            if ( c_i[2] != 0 && c_i[2] == - c_j[2] && c_i[0] == c_j[0] && c_i[1] == c_j[1] )
            {
                i_zinv[i] = j;
            }
        }
    }

    //--------------------------------------------------------------------------------------------//

    for (int y = 0; y < ny; y++ )
    {
        for (int x = 0; x < nx; x++ )
        {
            int *meio = ini_meio + x + y * nx + z * nx * ny;

            if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                for ( int i = 1; i < nvel; i++ )
                {
                    if ( i_zinv[i] != 0 )
                    {
                        double cup = f[i];

                        f[i]  = f[ i_zinv[i] ];

                        f[ i_zinv[i] ] = cup;
                    }
                }
            }
        }
    }
    
    delete[] i_zinv;
}

//================================================================================================//





//===================== Inverte a direção x (em um plano yz) =====================================//
//
//      Input: distribution function (*ini_f), geometry (*ini_meio), lattice vectors (*ini_c)
//              position of the plane xz, dimensions
//      Output:
//
//================================================================================================//

void inv_x ( double *ini_f, int *ini_meio, double *ini_c, int x, int nx, int ny, int nz )
{

//-------------- Determina as direções com x invertido ---------------------------------------//

    int *i_xinv = new int[nvel];

    for ( int i = 1; i < nvel; i++ )
    {
        i_xinv[i] = 0;

        for ( int j = i; j < nvel; j++ )
        {
            double *c_i = ini_c + i * dim;

            double *c_j = ini_c + j * dim;

            if ( c_i[0] != 0 && c_i[0] == - c_j[0] && c_i[2] == c_j[2] && c_i[1] == c_j[1] )
            {
                i_xinv[i] = j;
            }
        }
    }

    //--------------------------------------------------------------------------------------------//

    for (int y = 0; y < ny; y++ )
    {
        for (int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * nx * ny;

            if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                for ( int i = 1; i < nvel; i++ )
                {
                    if ( i_xinv[i] != 0 )
                    {
                        double cup = f[i];

                        f[i]  = f[ i_xinv[i] ];

                        f[ i_xinv[i] ] = cup;
                    }
                }
            }
        }
    }
    
    delete[] i_xinv;
}

//================================================================================================//




//===================== Inverte a função distribuição f[i] -> f[-i] ==============================//
//
//      Input: f[i]
//      Output:
//
//================================================================================================//

inline void inv_N ( double *N )
{
            
	for ( int i = 1; i < nvel; i = i + 2 )
	{
		double swap = N[i];
		
		N[i] = N[i+1];
		
		N[i+1] = swap;
	}

}
//================================================================================================//




//===================== Espelha o plano yz - direção x ===========================================//
//
//      Input: distribution function (*ini_f), geometry (*ini_meio), lattice vectors (*ini_c)
//              position of the plane xz, dimensions
//      Output:
//
//================================================================================================//

void bond_mirror_x ( double *ini_f, int *ini_meio, double *ini_c, int x, int nx, int ny, int nz )
{
    int infx;

    inv_x ( ini_f, ini_meio, ini_c, x, nx, ny, nz );

    if ( x < nx / 2 ) infx = 1;

    else infx = -1;

    int x_adj = x + infx;

    //-------------- Determina as direções a serem transferidas ----------------------------------//

    int *i_x = new int[nvel];

    for ( int i = 1; i < nvel; i++ )
    {
        double *c = ini_c + i * dim;

        if ( c[0] == infx )
        {
            i_x[i] = i;
        }
        else
        {
            i_x[i] = 0;
        }
    }

    //-------------- Transfere as distribuições para pontos adjacentes ---------------------------//

    for (int y = 0; y < ny; y++ )
    {
        for (int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * nx * ny;

            if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                int *meio_adj = ini_meio + x_adj + y * nx + z * nx * ny;

                if ( *meio_adj )
                {
                    double *f_adj = ini_f + ( *meio_adj - 1 ) * nvel;

                    for ( int i = 1; i < nvel; i++ )
                    {
                        if ( i_x[i] != 0 ) f_adj[ i_x[i] ] = f[ i_x[i] ];
                    }
                }
            }
        }
    }
    
    delete[] i_x;
}

//================================================================================================//




//===================== Espelha o plano xz - direção y ===========================================//
//
//      Input: distribution function (*ini_f), geometry (*ini_meio), lattice vectors (*ini_c)
//              position of the plane xz, dimensions
//      Output:
//
//================================================================================================//

void bond_mirror_y ( double *ini_f, int *ini_meio, double *ini_c, int y, int nx, int ny, int nz )
{
    int infy;

    inv_y ( ini_f, ini_meio, ini_c, y, nx, ny, nz );

    if ( y < ny / 2 ) infy = 1;

    else infy = -1;

    int y_adj = y + infy;

    //-------------- Determina as direções a serem transferidas ----------------------------------//

    int *i_y = new int[nvel];

    for ( int i = 1; i < nvel; i++ )
    {
        double *c = ini_c + i * dim;

        if ( c[1] == infy )
        {
            i_y[i] = i;
        }
        else
        {
            i_y[i] = 0;
        }
    }

    //-------------- Transfere as distribuições para pontos adjacentes ---------------------------//

    for (int z = 0; z < nz; z++ )
    {
        for (int x = 0; x < nx; x++ )
        {
            int *meio = ini_meio + x + y * nx + z * nx * ny;

            if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                int *meio_adj = ini_meio + x + y_adj * nx + z * nx * ny;

                if ( *meio_adj )
                {
                    double *f_adj = ini_f + ( *meio_adj - 1 ) * nvel;

                    for ( int i = 1; i < nvel; i++ )
                    {
                        if ( i_y[i] != 0 ) f_adj[ i_y[i] ] = f[ i_y[i] ];
                    }
                }
            }
        }
    }
    
    delete[] i_y;
}

//================================================================================================//




//===================== Espelha o plano xy - direção z ===========================================//
//
//      Input: distribution function (*ini_f), geometry (*ini_meio), lattice vectors (*ini_c)
//              position of the plane xz, dimensions
//      Output:
//
//================================================================================================//

void bond_mirror_z ( double *ini_f, int *ini_meio, double *ini_c, int z, int nx, int ny, int nz )
{
    int infz;

    inv_z ( ini_f, ini_meio, ini_c, z, nx, ny, nz );

    if ( z < nz / 2 ) infz = 1;

    else infz = -1;

    int z_adj = z + infz;

    //-------------- Determina as direções a serem transferidas ----------------------------------//

    int *i_z = new int[nvel];

    for ( int i = 1; i < nvel; i++ )
    {
        double *c = ini_c + i * dim;

        if ( c[2] == infz )
        {
            i_z[i] = i;
        }
        else
        {
            i_z[i] = 0;
        }
    }

    //-------------- Transfere as distribuições para pontos adjacentes ---------------------------//

    for (int y = 0; y < ny; y++ )
    {
        for (int x = 0; x < nx; x++ )
        {
            int *meio = ini_meio + x + y * nx + z * nx * ny;

            if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                int *meio_adj = ini_meio + x + y * nx + z_adj * nx * ny;

                if ( *meio_adj )
                {
                    double *f_adj = ini_f + ( *meio_adj - 1 ) * nvel;

                    for ( int i = 1; i < nvel; i++ )
                    {
                        if ( i_z[i] != 0 ) f_adj[ i_z[i] ] = f[ i_z[i] ];
                    }
                }
            }
        }
    }
    
    delete[] i_z;
}

//================================================================================================//




//===================== Acrescenta força em um sítio =============================================//
//
//      Input: distribution function, lattice vectors, force in the x,y,z directions
//      Output: force efectivelly imposed
//
//================================================================================================//

double force ( double *f, double *ini_c, double gx, double gy, double gz )
{

    double sum_force = 0.0;
    double g[dim];
    double *c;
    double gforce;

    g[0] = gx;
    g[1] = gy;
    g[2] = gz;

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        gforce = prod_int ( g, c );

        f[i] = f[i] + gforce;

        sum_force = sum_force + fabs ( gforce );
    }

    return sum_force;
}

//================================================================================================//




//===================== Grava os o meio em formato .vtk ==========================================//
//
//      Input: geometry, name, dimensions
//      Output:
//
//================================================================================================//

void rec_geometria ( int *ini_meio, char *nome, int nx, int ny, int nz )
{
    int *meio;

    ofstream fmeio_out( nome );

    fmeio_out << "# vtk DataFile Version 2.0" << endl;
    fmeio_out << "Geometria" << endl;
    fmeio_out << "ASCII" << endl;
    fmeio_out << "DATASET STRUCTURED_POINTS" << endl;
    fmeio_out << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fmeio_out << "ASPECT_RATIO 1 1 1" << endl;
    fmeio_out << "ORIGIN 0 0 0" << endl;
    fmeio_out << "POINT_DATA " << nx * ny * nz << endl;
    fmeio_out << "SCALARS geometria int" << endl;
    fmeio_out << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio ) fmeio_out << 1 << " ";

                else fmeio_out << 0 << " ";
            }
            fmeio_out << endl;
        }
    }

    fmeio_out.close();

}

//================================================================================================//




//===================== Grava os resultados - monofásico =========================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_one_fluid ( int *ini_meio, double *ini_f, double *ini_c, double rho_ini,
                     unsigned int passo, int nx, int ny, int nz )
{
    double rho, vx, vy, vz;
    double *f;
    int *meio;

    char nomero[50];
    char nomevel[50];

    sprintf ( nomero, "rho_%06d.vtk", passo );
    sprintf ( nomevel, "vel_%06d.vtk", passo );

    ofstream fro ( nomero );
    ofstream fvel ( nomevel );

    fro << "# vtk DataFile Version 2.0" << endl;
    fro << "Densidade" << endl;
    fro << "ASCII" << endl;
    fro << "DATASET STRUCTURED_POINTS" << endl;
    fro << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fro << "ASPECT_RATIO 1 1 1" << endl;
    fro << "ORIGIN 0 0 0" << endl;
    fro << "POINT_DATA " << nx*ny*nz << endl;
    fro << "SCALARS densidade double" << endl;
    fro << "LOOKUP_TABLE default" << endl;

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;

                    calcula ( f, ini_c, &vx, &vy, &vz, &rho );

                    fro << rho << " ";
                    fvel << vx  << " " << vy << " " << vz << " ";
                }
                else
                {
                    fro << 0.0 << " ";
                    fvel << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }

            fro << endl;
            fvel << endl;
        }

    }

    fro.close();
    fvel.close();

}

//================================================================================================//





//===================== Grava o campo de velocidade - monofásico =================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_velocity ( int *ini_meio, double *ini_f, double *ini_c, unsigned int passo, int nx,
                    int ny, int nz )
{
    char nomevel[50];

    sprintf ( nomevel, "vel_%06d.vtk", passo );

    ofstream fvel ( nomevel );

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;
    
    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int *meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    double *f = ini_f + ( *meio - 1 ) * nvel;

                    double vx, vy, vz, rho;

                    calcula ( f, ini_c, &vx, &vy, &vz, &rho );

                    fvel << vx  << " " << vy << " " << vz << " ";
                }
                else
                {
                    fvel << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }
            fvel << endl;
        }

    }
    fvel.close();
}

//================================================================================================//




//===================== Grava o campo de velocidade - monofásico =================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_velocity_bin ( int *ini_meio, double *ini_f, double *ini_c, unsigned int passo, int nx,
                    int ny, int nz )
{
    char nomevel[50];

    sprintf ( nomevel, "vel_%06d.vtk", passo );

    ofstream fvel;
    
    fvel.open ( nomevel );

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "BINARY" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "SPACING 1 1 1" << endl;    
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;
    
    fvel.close();
    
    fvel.open ( nomevel, ios::app | ios::binary );

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
				double vx = 0.0;
				double vy = 0.0;
				double vz = 0.0;
				
				double rho = 0.0;				
				
                int *meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    double *f = ini_f + ( *meio - 1 ) * nvel;

                    calcula ( f, ini_c, &vx, &vy, &vz, &rho );
                    
					fvel.write ( ( char* ) &vx, sizeof ( double ) );
					fvel.write ( ( char* ) &vx, sizeof ( double ) );
					fvel.write ( ( char* ) &vy, sizeof ( double ) );
					fvel.write ( ( char* ) &vz, sizeof ( double ) );
                }
                else
                {                   
                    fvel.write ( ( char* ) &vx, sizeof ( double ) );
					fvel.write ( ( char* ) &vy, sizeof ( double ) );
					fvel.write ( ( char* ) &vz, sizeof ( double ) );
                }
            }
            fvel << endl;
        }

    }
    fvel.close();
}

//================================================================================================//



//===================== Grava o campo de velocidade - dois fluidos ===============================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_velocity_twofluid ( int *ini_meio, double *ini_R, double *ini_B, double *ini_c, int passo,
                             int nx, int ny, int nz )
{
    char nomevel[50];

    sprintf ( nomevel, "vel_%06d.vtk", passo );

    ofstream fvel ( nomevel );

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int *meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    double *R = ini_R + ( *meio - 1 ) * nvel;

                    double *B = ini_B + ( *meio - 1 ) * nvel;

                    double vxR, vyR, vzR, rhoR;

                    calcula ( R, ini_c, &vxR, &vyR, &vzR, &rhoR );

                    double vxB, vyB, vzB, rhoB;

                    calcula ( B, ini_c, &vxB, &vyB, &vzB, &rhoB );

                    double concR = rhoR / ( rhoR + rhoB );

                    double concB = 1.0 - concR;

                    double vx = concR * vxR + concB * vxB;

                    double vy = concR * vyR + concB * vyB;

                    double vz = concR * vzR + concB * vzB;

                    fvel << vx  << " " << vy << " " << vz << " ";
                }
                else
                {
                    fvel << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }
            fvel << endl;
        }

    }
    fvel.close();
}

//================================================================================================//




//===================== Grava o campo de densidade - um fluido ===================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_density ( int *ini_meio, double *ini_f, double *ini_c, unsigned int passo, int nx,
                   int ny, int nz, string nomerho )
{
    nomerho = nomerho + to_string( passo ) + ".vtk";

    ofstream frho ( nomerho );

    frho << "# vtk DataFile Version 2.0" << endl;
    frho << "Densidade" << endl;
    frho << "ASCII" << endl;
    frho << "DATASET STRUCTURED_POINTS" << endl;
    frho << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    frho << "ASPECT_RATIO 1 1 1" << endl;
    frho << "ORIGIN 0 0 0" << endl;
    frho << "POINT_DATA " << nx*ny*nz << endl;
    frho << "SCALARS densidade double" << endl;
    frho << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int *meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    double *f = ini_f + ( *meio - 1 ) * nvel;

                    double rho = mass( f );

                    frho << rho << " ";
                }
                else
                {
                    frho << 0.0 << " ";
                }
            }

            frho << endl;
        }

    }
    frho.close();
}
//================================================================================================//




//===================== Grava os resultados - monofásico - térmico ===============================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_one_fluid_thermal ( int *ini_meio, double *ini_f, double *ini_c, double rho_ini,
                             unsigned int passo, int nx, int ny, int nz )
{
    double rho, vx, vy, vz, theta;
    double *f;
    int *meio;

    char nomero[50];
    char nomevel[50];

    sprintf ( nomero, "rho_%06d.vtk", passo );
    sprintf ( nomevel, "vel_%06d.vtk", passo );

    ofstream fro ( nomero );
    ofstream fvel ( nomevel );

    fro << "# vtk DataFile Version 2.0" << endl;
    fro << "Densidade" << endl;
    fro << "ASCII" << endl;
    fro << "DATASET STRUCTURED_POINTS" << endl;
    fro << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fro << "ASPECT_RATIO 1 1 1" << endl;
    fro << "ORIGIN 0 0 0" << endl;
    fro << "POINT_DATA " << nx*ny*nz << endl;
    fro << "SCALARS densidade double" << endl;
    fro << "LOOKUP_TABLE default" << endl;

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;

                    calcula_thermal ( f, ini_c, &vx, &vy, &vz, &rho, &theta );

                    fro << rho << " ";
                    fvel << vx  << " " << vy << " " << vz << " ";
                }
                else
                {
                    fro << 0.0 << " ";
                    fvel << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }

            fro << endl;
            fvel << endl;
        }

    }

    fro.close();
    fvel.close();

}

//================================================================================================//




//===================== Calcula a vorticidade em um ponto x, y, z ================================//
//
//      Input: geometry, distributions functions, lattice vectors, position, dimensions
//      Output: vorticity
//
//================================================================================================//

void vorticity ( int *ini_meio, double *ini_c, double *ini_f, int x, int y, int z,
                 int nx, int ny, int nz, double *vort_x, double *vort_y, double *vort_z )

{
    double *f, rho;

    int x_p = ( x + 1 ) % nx;
    int x_m = ( x - 1 + nx ) % nx;

    int y_p = ( y + 1 ) % ny;
    int y_m = ( y - 1 + ny ) % ny;

    int z_p = ( z + 1 ) % nz;
    int z_m = ( z - 1 + nz ) % nz;

    int *x_mais, *x_menos, *y_mais, *y_menos, *z_mais, *z_menos;

    x_mais = ini_meio + x_p + y * nx + z * ny * nx;
    x_menos = ini_meio + x_m + y * nx + z * ny * nx;

    y_mais = ini_meio + x + y_p * nx + z * ny * nx;
    y_menos = ini_meio + x + y_m * nx + z * ny * nx;

    z_mais = ini_meio + x + y * nx + z_p * ny * nx;
    z_menos = ini_meio + x + y * nx + z_m * ny * nx;

    double vx_x_mais, vx_x_menos, vy_x_mais, vy_x_menos, vz_x_mais, vz_x_menos;

    if ( *x_mais )
    {
        f = ini_f + ( *x_mais - 1 ) * nvel;

        calcula ( f, ini_c, &vx_x_mais, &vy_x_mais, &vz_x_mais, &rho );
    }
    else
    {
        vx_x_mais = vy_x_mais = vz_x_mais = 0.0;
    }

    if ( *x_menos )
    {
        f = ini_f + ( *x_menos - 1 ) * nvel;

        calcula ( f, ini_c, &vx_x_menos, &vy_x_menos, &vz_x_menos, &rho );
    }
    else
    {
        vx_x_menos = vy_x_menos = vz_x_menos = 0.0;
    }

    double vx_y_mais, vx_y_menos, vy_y_mais, vy_y_menos, vz_y_mais, vz_y_menos;

    if ( *y_mais )
    {
        f = ini_f + ( *y_mais - 1 ) * nvel;

        calcula ( f, ini_c, &vx_y_mais, &vy_y_mais, &vz_y_mais, &rho );
    }
    else
    {
        vx_y_mais = vy_y_mais = vz_y_mais = 0.0;
    }

    if ( *y_menos )
    {
        f = ini_f + ( *y_menos - 1 ) * nvel;

        calcula ( f, ini_c, &vx_y_menos, &vy_y_menos, &vz_y_menos, &rho );
    }
    else
    {
        vx_y_menos = vy_y_menos = vz_y_menos = 0.0;
    }

    double vx_z_mais, vx_z_menos, vy_z_mais, vy_z_menos, vz_z_mais, vz_z_menos;

    if ( *z_mais )
    {
        f = ini_f + ( *z_mais - 1 ) * nvel;

        calcula ( f, ini_c, &vx_z_mais, &vy_z_mais, &vz_z_mais, &rho );
    }
    else
    {
        vx_z_mais = vy_z_mais = vz_z_mais = 0.0;
    }

    if ( *z_menos )
    {
        f = ini_f + ( *z_menos - 1 ) * nvel;

        calcula ( f, ini_c, &vx_z_menos, &vy_z_menos, &vz_z_menos, &rho );
    }
    else
    {
        vx_z_menos = vy_z_menos = vz_z_menos = 0.0;
    }

    double deriv_vz_y = ( vz_y_mais - vz_y_menos ) / 2.0;
    double deriv_vy_z = ( vy_z_mais - vy_z_menos ) / 2.0;

    *vort_x = deriv_vz_y - deriv_vy_z;

    double deriv_vx_z = ( vx_z_mais - vx_z_menos ) / 2.0;
    double deriv_vz_x = ( vz_x_mais - vz_x_menos ) / 2.0;

    *vort_y = deriv_vx_z - deriv_vz_x;

    double deriv_vy_x = ( vy_x_mais - vy_x_menos ) / 2.0;
    double deriv_vx_y = ( vx_y_mais - vx_y_menos ) / 2.0;

    *vort_z = deriv_vy_x - deriv_vx_y;
}

//================================================================================================//





//===================== Grava os resultados - vorticidade ========================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_vorticity ( int *ini_meio, double *ini_f, double *ini_c, unsigned int passo, int nx,
                     int ny, int nz )
{
    int *meio;

    char nome_vorticity[50];

    sprintf ( nome_vorticity, "vor_%06d.vtk", passo );

    ofstream fvor ( nome_vorticity );

    fvor << "# vtk DataFile Version 2.0" << endl;
    fvor << "Velocidade" << endl;
    fvor << "ASCII" << endl;
    fvor << "DATASET STRUCTURED_POINTS" << endl;
    fvor << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvor << "ASPECT_RATIO 1 1 1" << endl;
    fvor << "ORIGIN 0 0 0" << endl;
    fvor << "POINT_DATA " << nx*ny*nz << endl;
    fvor << "VECTORS vorticidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y * nx + z * ny * nx;

                if ( *meio )
                {
                    double vort_x, vort_y, vort_z;

                    vorticity ( ini_meio, ini_c, ini_f, x, y, z, nx, ny, nz,
                                &vort_x, &vort_y, &vort_z );

                    fvor << vort_x  << " " << vort_y << " " << vort_z << " ";
                }
                else
                {
                    fvor << 0.0  << " " << 0.0 << " " << 0.0 << " ";
                }
            }

            fvor << endl;
        }

    }

    fvor.close();
}

//================================================================================================//




//===================== Grava arquivo de recuperação de dados ====================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_recovery ( double *ini_f, double *ini_c, int ptos_meio, int passo, int precision )
{
    cout << "\nGravando arquivo de recuperação..." << endl;

    ofstream f_rec ( "Arq_rec.dat", ios::binary );

    f_rec << passo << endl;

    for ( int pto = 0; pto < ptos_meio; pto++ )
    {
        double *f = ini_f + ( pto ) * nvel;

        double vx, vy, vz, rho;

        calcula ( f, ini_c, &vx, &vy, &vz, &rho );

        f_rec << setprecision( precision ) << vx << " ";

        f_rec << setprecision( precision ) << vy << " ";

        f_rec << setprecision( precision ) << vz << " ";

        f_rec << setprecision( precision ) << rho << " ";

    }

    f_rec.close();

    cout << "... ... ... !" << endl << endl;
}

//================================================================================================//




//===================== Grava arquivo de recuperação de dados ====================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_recovery_two_fluids ( double *ini_f, double *ini_g, double *ini_c, int ptos_meio, 
								int passo )
{
	long int tamanho = nvel * ptos_meio * sizeof ( double );
	
    cout << "\nGravando arquivo de recuperação..." << endl;

    fstream f_rec ( "Arq_rec.dat", ios::out | ios::binary );
    
	f_rec.write ( ( char* ) &passo, sizeof (int) );	
	
	f_rec.write ( ( char* ) ini_f, tamanho );
	
	f_rec.write ( ( char* ) ini_g, tamanho );

	f_rec.close();
	
	cout << "... ... ... !" << endl << endl;
}

//================================================================================================//





//===================== Condição de contorno de derivada nula (rho e vel.) - direção x ===========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_x ( int *ini_meio , double *ini_f, double *ini_c, int posx,
                       int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y*nx + z*ny*nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA;
                
                double rhoA = 0.0;

                double *f;

                meio = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y*nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vxA, vyA, vzA, rhoA, ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
//================================================================================================//




//===================== Condição de contorno de derivada nula da velocidade - direção x ==========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//
/*/
void bond_eqdvnull_rho_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double rho,
                           int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA, rhoA;

                double *f;

                meio = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y*nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vxA, vyA, vzA, rho, ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
/*/
//================================================================================================//




//===================== Condição de contorno de derivada nula da velocidade - direção x ==========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_dvnull_rho_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double rho,
                           int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int x = posx;
    
    int infx = -1;
    
    if ( posx < nx / 2 ) infx = 1;

	#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int* meio_pto = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = ini_f + ( *meio_pto - 1 ) * nvel;
				
                int* meio_adj = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                if ( *meio_adj )
                {
                    double* f_adj = ini_f + ( *meio_adj - 1 ) * nvel;
                    
                    double rho_adj = mass( f_adj );
                    
                    double fator = rho / rho_adj;
                
					for ( int i = 0; i < nvel; i++ ) f_pto[i] = f_adj[i] * fator;
                }
                else
                {
                    dist_eq ( 0., 0., 0., rho, ini_c, f_pto, W, one_over_c_s2 );
                }
            }
        }
    }
}
//================================================================================================//



//===================== Condição de contorno de derivada nula - direção x ========================//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_dvnull_x_d3q19 ( int *ini_meio , double *ini_f, double *ini_c, int posx,
                           int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int x = posx;
    
    int infx = -1;
    
    if ( posx < nx / 2 ) infx = 1;

	#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int* meio_pto = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio_pto )
            {
				double* f_pto = ini_f + ( *meio_pto - 1 ) * nvel;
				
                int* meio_adj = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                if ( *meio_adj && infx < 0 )
                {
                    double* f_adj = ini_f + ( *meio_adj - 1 ) * nvel; 
                                    
					f_pto[2] = f_adj[2];
					f_pto[8] = f_adj[8];
					f_pto[10] = f_adj[10];
					f_pto[12] = f_adj[12];
					f_pto[14] = f_adj[14];
					
                }
                else if ( *meio_adj && infx > 0 )
                {
                    double* f_adj = ini_f + ( *meio_adj - 1 ) * nvel; 
                                    
					f_pto[1] = f_adj[1];
					f_pto[7] = f_adj[7];
					f_pto[9] = f_adj[9];
					f_pto[11] = f_adj[11];
					f_pto[13] = f_adj[13];
					
                }
                else
                {
                    //dist_eq ( 0., 0., 0., rho, ini_c, f_pto, W, one_over_c_s2 );
                }
            }
        }
    }
}
//================================================================================================//




//===================== Condição de contorno de derivada nula da velocidade - direção x ==========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_cbc_x ( int *ini_meio, double *ini_f, double *ini_f_old, double *ini_c, int posx, 
                           int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula a velocidade -----------------------------------------------//
                
				int* meio_A = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                double *f_A = ini_f + ( *meio_A - 1 ) * nvel;
                
                double U, vy, vz, rho;
                
                calcula ( f_A, ini_c, &U, &vy, &vz, &rho );
                
                //--------------------------------------------------------------------------------//
                                
                meio = ini_meio + x + y * nx + z * ny * nx;

                double *f = ini_f + ( *meio - 1 ) * nvel;
                
                double *f_old = ini_f_old + ( *meio - 1 ) * nvel;
                
                double factor = 1.0 / ( 1.0 + U );
                
                for ( int i = 0 ; i < nvel; i++ )
				{
					f[i] = factor * ( f_old[i] + U * f_A[i]);
				}

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
//================================================================================================//




//===================== Condição de contorno de convectiva (bifásica) - direção x ================//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_cbc_2fase_x ( int *ini_meio, double *ini_f, double *ini_f_old, double *ini_f_other, 
						double *ini_c, int posx, int nx, int ny, int nz, 
						double *W, double one_over_c_s2 )
{
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula a velocidade -----------------------------------------------//
                
				int* meio_A = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                double *f_A = ini_f + ( *meio_A - 1 ) * nvel;
                
                double U_one, vy, vz, rho_one;
                
                calcula ( f_A, ini_c, &U_one, &vy, &vz, &rho_one );
                
                double *f_A_other = ini_f_other + ( *meio_A - 1 ) * nvel;
                
                double U_two, rho_two;
                
                calcula ( f_A_other, ini_c, &U_two, &vy, &vz, &rho_two );
                
                double U = ( U_one * rho_one + U_two * rho_two ) / ( rho_one + rho_two );
                
                //--------------------------------------------------------------------------------//
                                
                meio = ini_meio + x + y * nx + z * ny * nx;

                double *f = ini_f + ( *meio - 1 ) * nvel;
                
                double *f_old = ini_f_old + ( *meio - 1 ) * nvel;
                
                double factor = 1.0 / ( 1.0 + U );
                
                for ( int i = 0 ; i < nvel; i++ )
				{
					f[i] = factor * ( f_old[i] + U * f_A[i]);
				}

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
//================================================================================================//





//===================== Condição de contorno de convectiva (mediadores) - direção x ==============//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_cbc_med_x ( int *ini_meio, double *ini_f, double *ini_f_old, double *ini_R, 
						double *ini_B, double *ini_c, int posx, int nx, int ny, int nz, 
						double *W, double one_over_c_s2 )
{
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula a velocidade -----------------------------------------------//
                
				int* meio_A = ini_meio + ( x + infx ) + y * nx + z * ny * nx;
				
				double *f_A = ini_f + ( *meio_A - 1 ) * nvel;

                double *R_A = ini_R + ( *meio_A - 1 ) * nvel;
                
                double U_red, vy, vz, rho_red;
                
                calcula ( R_A, ini_c, &U_red, &vy, &vz, &rho_red );
                
                double *B_A = ini_B + ( *meio_A - 1 ) * nvel;
                
                double U_blue, rho_blue;
                
                calcula ( B_A, ini_c, &U_blue, &vy, &vz, &rho_blue );
                
                double U = ( U_red * rho_red + U_blue * rho_blue ) / ( rho_red + rho_blue );
                
                //--------------------------------------------------------------------------------//
                                
                meio = ini_meio + x + y * nx + z * ny * nx;

                double *f = ini_f + ( *meio - 1 ) * nvel;
                
                double *f_old = ini_f_old + ( *meio - 1 ) * nvel;
                
                double factor = 1.0 / ( 1.0 + U );
                
                for ( int i = 0 ; i < nvel; i++ )
				{
					f[i] = factor * ( f_old[i] + U * f_A[i]);
				}

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
//================================================================================================//





//===================== Condição de contorno de derivada nula da velocidade - direção x ==========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_rholocal_x ( int *ini_meio , double *ini_f, double *ini_c, int posx,
                                int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y*nx + z*ny*nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA, rhoA;

                double *f;

                meio = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y*nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;


                dist_eq ( vxA, vyA, vzA, mass( f ), ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
//================================================================================================//



//===================== Condição de contorno de derivada nula da velocidade - direção x ==========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_rholocalveldim_x ( int *ini_meio , double *ini_f, double *ini_c, int posx,
                                      double fat_at, int nx, int ny, int nz, double *W, 
                                      double one_over_c_s2 )
{
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y*nx + z*ny*nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA, rhoA;

                double *f;

                meio = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }

                //------------- Atenuação da velocidade ------------------------------------------//

                double vx_imp = fat_at * vxA;

                double vy_imp = fat_at * vyA;

                double vz_imp = fat_at * vzA;

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y*nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vx_imp, vy_imp, vz_imp, mass( f ), ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
//================================================================================================//



//===================== Condição de contorno de derivada nula da velocidade - direção x ==========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eq_rhovel_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double rho,
                        double vx, double vy, double vz, int nx,  int ny, int nz, double *W,
                        double one_over_c_s2  )
{

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + posx + y*nx + z * ny * nx;

            if ( *meio )
            {
                //--------------------------------------------------------------------------------//

                double *f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vx, vy, vz, rho, ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno Zou & He - direção x ================================//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_ZouHeD3Q19_rhovel_x ( int *ini_meio , double *ini_f, int posx, double rho,
                                double vx, double vy, double vz, int nx, int ny, int nz, double *W,
                                double one_over_c_s2 )
{
    int infx = 0;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio && infx )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                f[1]  = f[1] + rho * vx / 3.0;

                f[7]  = f[7] + rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[9]  = f[9] + rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[11] = f[11] + rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[13] = f[13] + rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
            else if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                f[2]  = f[2] - rho * vx / 3.0;

                f[10] = f[10] - rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[8]  = f[8] - rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[14] = f[14] - rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[12] = f[12] - rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
        }
    }
}
//================================================================================================//




//===================== Condição de contorno Zou & He - rho imposto - direção x ==================//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_ZouHeD3Q19_vx_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double vx,
                            int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int x = posx;

    int infx = 0;

    if ( posx < nx / 2 ) infx = 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio && infx )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                //------------------- Calcula a soma das incógnitas ------------------------------//

                double sum_f_xminus = f[1] + f[7] + f[9] + f[11] + f[13];

                double sum_f_xzero = f[0] + f[3] + f[4] + f[5] + f[6]
                                     + f[15] + f[16] + f[17] + f[18];

                //-------------------- Determinação de vx e imposição de vy e vz -----------------//

                double rho = ( 2.0 * sum_f_xminus + sum_f_xzero ) / ( 1.0 - vx );

                double vy = 0.0;

                double vz = 0.0;

                //--------------------------------------------------------------------------------//

                f[1]  = f[1] + rho * vx / 3.0;

                f[7]  = f[7] + rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[9]  = f[9] + rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[11] = f[11] + rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[13] = f[13] + rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
            else if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                //------------------- Calcula a soma das incógnitas ------------------------------//

                double sum_f_xplus = f[2] + f[8] + f[10] + f[12] + f[14];

                double sum_f_xzero = f[0] + f[3] + f[4] + f[5] + f[6]
                                     + f[15] + f[16] + f[17] + f[18];

                //-------------------- Determinação de vx e imposição de vy e vz -----------------//

                double rho = ( 2.0 * sum_f_xplus + sum_f_xzero ) / ( 1.0 - vx );

                double vy = 0.0;

                double vz = 0.0;

                //--------------------------------------------------------------------------------//

                f[2]  = f[2] - rho * vx / 3.0;

                f[10] = f[10] - rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[8]  = f[8] - rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[14] = f[14] - rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[12] = f[12] - rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno Zou & He - vel. local direção x =====================//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_ZouHeD3Q19_rho_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double rho,
                             int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int x = posx;

    int infx = 0;

    if ( posx < nx / 2 ) infx = 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio && infx )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                //------------------- Calcula a soma das incógnitas ------------------------------//

                double sum_f_xminus = f[1] + f[7] + f[9] + f[11] + f[13];

                double sum_f_xzero = f[0] + f[3] + f[4] + f[5] + f[6]
                                     + f[15] + f[16] + f[17] + f[18];

                double sum_f_xplus = rho - sum_f_xminus - sum_f_xzero;

                //-------------------- Determinação de vx e imposição de vy e vz -----------------//

                double vx = ( sum_f_xplus - sum_f_xminus ) / ( rho );

                double vy = 0.0;

                double vz = 0.0;

                //--------------------------------------------------------------------------------//

                f[1]  = f[1] + rho * vx / 3.0;

                f[7]  = f[7] + rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[9]  = f[9] + rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[11] = f[11] + rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[13] = f[13] + rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
            else if ( *meio )
            {
                double *f = ini_f + ( *meio - 1 ) * nvel;

                //------------------- Calcula a soma das incógnitas ------------------------------//

                double sum_f_xplus = f[2] + f[8] + f[10] + f[12] + f[14];

                double sum_f_xzero = f[0] + f[3] + f[4] + f[5] + f[6]
                                     + f[15] + f[16] + f[17] + f[18];

                double sum_f_xminus = rho - sum_f_xplus - sum_f_xzero;

                //-------------------- Determinação de vx e imposição de vy e vz -----------------//

                double vx = ( sum_f_xplus - sum_f_xminus ) / rho;

                double vy = 0.0;

                double vz = 0.0;

                //--------------------------------------------------------------------------------//

                f[2]  = f[2] - rho * vx / 3.0;

                f[10] = f[10] - rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[8]  = f[8] - rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[14] = f[14] - rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[12] = f[12] - rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
        }
    }
}

//================================================================================================//




//===================== Cond. de contorno de equilíbrio realimentado (rho & vel.) - x ============//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, velocities,
//             dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqretro_rhovel_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double rho,
                             double vx, double vy, double vz, int nx, int ny, int nz, double *W,
                             double one_over_c_s2 )
{

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + posx + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula densidade e velocidades do ponto ---------------------------//

                double vxP, vyP, vzP, rhoP;

                double *f = ini_f + ( *meio - 1 ) * nvel;

                calcula ( f, ini_c, &vxP, &vyP, &vzP, &rhoP );

                //--------------------------------------------------------------------------------//

                double rho_imp = rho - 0.02 * ( rhoP - rho ) ;

                double vx_imp = vx - 0.02 * ( vxP - vx );

                double vy_imp = vy - 0.02 * ( vyP - vy );

                double vz_imp = vz - 0.02 * ( vzP - vz );

                //--------------------------------------------------------------------------------//

                dist_eq ( vx_imp, vy_imp, vz_imp, rho_imp, ini_c, f, W, one_over_c_s2 );
            }
        }
    }
}

//================================================================================================//




//===================== Cond. de contorno de equilíbrio realimentado (rho & vel.) - y  ===========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, velocities,
//             dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqretro_rhovel_y ( int *ini_meio , double *ini_f, double *ini_c, int posy, double rho,
                             double vx, double vy, double vz, int nx, int ny, int nz, double *W,
                             double one_over_c_s2 )
{

#pragma omp parallel for

    for ( int x = 0; x < nx; x++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + posy * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula densidade e velocidades do ponto ---------------------------//

                double vxP, vyP, vzP, rhoP;

                double *f = ini_f + ( *meio - 1 ) * nvel;

                calcula ( f, ini_c, &vxP, &vyP, &vzP, &rhoP );

                //--------------------------------------------------------------------------------//

                double rho_imp = rho - 0.02 * ( rhoP - rho ) ;

                double vx_imp = vx - 0.02 * ( vxP - vx );

                double vy_imp = vy - 0.02 * ( vyP - vy );

                double vz_imp = vz - 0.02 * ( vzP - vz );

                //--------------------------------------------------------------------------------//

                dist_eq ( vx_imp, vy_imp, vz_imp, rho_imp, ini_c, f, W, one_over_c_s2 );
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno de derivada nula (Neumann) - direração y ============//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_rho_y ( int *ini_meio , double *ini_f, double *ini_c, int posy, double rho, 
							int nx, int ny, int nz, double *W, double one_over_c_s2 )
{

    int infy;
    int y = posy;

    if ( posy < ny / 2 ) infy = 1;
    else infy = - 1;

#pragma omp parallel for

    for ( int x = 0; x < nx; x++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA, rhoA;

                double *f;

                meio = ini_meio + x + ( y + infy ) * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y * nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vxA, vyA, vzA, rho, ini_c, f, W, one_over_c_s2 );
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno de derivada nula da densidade - direção x ===========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_vx_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double vx,
                           int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int infx;
    int x = posx;

    if ( posx < nx / 2 ) infx = 1;
    else infx = - 1;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA, rhoA;

                double *f;

                meio = ini_meio + ( x + infx ) + y * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    rhoA = vxA = vyA = vzA = 0.0;
                }

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y*nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vx, vyA, vzA, rhoA, ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
//================================================================================================//





//===================== Condição de contorno de derivada nula (Neumann) - direração y ============//
//
//      Input: geometry, distribution functions, lattice vectors, position, vx,vy,vz, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_vy_y ( int *ini_meio , double *ini_f, double *ini_c, int posy, double vy, int nx,
                           int ny, int nz, double *W, double one_over_c_s2 )
{

    int infy;
    int y = posy;

    if ( posy < ny / 2 ) infy = 1;
    else infy = - 1;

#pragma omp parallel for

    for ( int x = 0; x < nx; x++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA;
                double rhoA = 0.0;

                double *f;

                meio = ini_meio + x + ( y + infy ) * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y * nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vxA, vy, vzA, rhoA, ini_c, f, W, one_over_c_s2 );
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno de derivada nula (Neumann) - direração y ============//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_rholocal_y ( int *ini_meio , double *ini_f, double *ini_c, int posy,
                                int nx, int ny, int nz, double *W, double one_over_c_s2 )
{

    int infy;
    int y = posy;

    if ( posy < ny / 2 ) infy = 1;
    else infy = - 1;

#pragma omp parallel for

    for ( int x = 0; x < nx; x++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA, rhoA;

                double *f;

                meio = ini_meio + x + ( y + infy ) * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y * nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vxA, vyA, vzA, mass( f ), ini_c, f, W, one_over_c_s2 );
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno de derivada nula (Neumann) - direração y ============//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_rholocalveldim_y ( int *ini_meio , double *ini_f, double *ini_c, int posy,
                                      double fat_at, int nx, int ny, int nz, double *W, 
                                      double one_over_c_s2 )
{

    int infy;
    int y = posy;

    if ( posy < ny / 2 ) infy = 1;
    else infy = - 1;

#pragma omp parallel for

    for ( int x = 0; x < nx; x++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA, rhoA;

                double *f;

                meio = ini_meio + x + ( y + infy ) * nx + z * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;
                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }

                //------------- Atenuação da velocidade ------------------------------------------//

                double vx_imp = fat_at * vxA;

                double vy_imp = fat_at * vyA;

                double vz_imp = fat_at * vzA;

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y * nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vx_imp, vy_imp, vz_imp, mass( f ), ini_c, f, W, one_over_c_s2 );
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno de equilíbrio - direção y ===========================//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, velocities,
//             dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eq_rhovel_y ( int *ini_meio , double *ini_f, double *ini_c, int posy, double rho, 
						double vx, double vy, double vz, int nx,  int ny, int nz, double *W, 
                        double one_over_c_s2 )
{
    int y = posy;

#pragma omp parallel for

    for ( int x = 0; x < nx; x++ )
    {
        for ( int z = 0; z < nz; z++ )
        {
            int *meio = ini_meio + x + y*nx + z*ny*nx;

            if ( *meio )
            {
                //--------------------------------------------------------------------------------//

                double *f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vx, vy, vz, rho, ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno de derivada nula (Neumann) - direração z ============//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqdvnull_rho_z ( int *ini_meio, double *ini_f,  double *ini_c, int posz, double rho,
                           int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    int infz;
    int z = posz;

    if ( posz < nz / 2 ) infz = 1;
    else infz = - 1;

#pragma omp parallel for

    for ( int x = 0; x < nx; x++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //----------- Calcula velocidade do ponto adjacente ------------------------------//

                double vxA, vyA, vzA, rhoA;

                double *f;

                meio = ini_meio + x + y * nx + ( z + infz ) * ny * nx;

                if ( *meio )
                {
                    f = ini_f + ( *meio - 1 ) * nvel;

                    calcula ( f, ini_c, &vxA, &vyA, &vzA, &rhoA );
                }
                else
                {
                    vxA = vyA = vzA = 0.0;
                }  

                //--------------------------------------------------------------------------------//

                meio = ini_meio + x + y * nx + z * ny * nx;

                f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vxA, vyA, vzA, rho, ini_c, f, W, one_over_c_s2 );
            }
        }
    }
}

//================================================================================================//




//===================== Condição de contorno de equilíbrio - direção z ===========================//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, velocities,
//             dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eq_rhovel_z ( int *ini_meio , double *ini_f, double *ini_c, int posz, 
						double rho, double vx, double vy, double vz, int nx,  int ny, int nz, 
						double *W, double one_over_c_s2 )
{
    int z = posz;

#pragma omp parallel for

    for ( int y = 0; y < ny; y++ )
    {
        for ( int x = 0; x < nx; x++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //--------------------------------------------------------------------------------//

                double *f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vx, vy, vz, rho, ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }
        }
    }
}
//================================================================================================//




//===================== Condição de contorno de equilíbrio - fluxo parabólico - direção x ========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, velocities,
//             dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_eqparabolly_x ( int *ini_meio , double *ini_f, double *ini_c, int posx, double rho,
                          double vx_max, int L, int nx, int ny, int nz, double *W,
                          double one_over_c_s2 )
{
    int x = posx;

    double vy = 0.0;

    double vz = 0.0;

    L++;

    double *vel = new double[ L ];

    for ( int h = 0; h < L; h++ )
    {
        double d_h = ( double ) h;

        vel[ h ] = -4.0 * vx_max / ( L * L ) * ( d_h ) * ( d_h - L );
    }

#pragma omp parallel for

    for ( int z = 0; z < nz; z++ )
    {
        int h = 0;

        for ( int y = 0; y < ny; y++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //--------------------------------------------------------------------------------//

                h++;

                double vx = vel[ h ];

                double *f = ini_f + ( *meio - 1 ) * nvel;

                dist_eq ( vx, vy, vz, rho, ini_c, f, W, one_over_c_s2 );

                //--------------------------------------------------------------------------------//
            }

        }
    }

    delete[] vel;
}
//================================================================================================//




//===================== Condição de contorno de equilíbrio - fluxo parabólico - direção x ========//
//
//      Input: geometry, distribution functions, lattice vectors, position, density, velocities,
//             dimensions
//      Output: distribution function
//
//================================================================================================//

void bond_parabolly_x ( int *ini_meio , double *ini_f, double *ini_c, int posx,
                          double vx_max, int L, int nx, int ny, int nz, double *W,
                          double one_over_c_s2 )
                       
{
    int x = posx;
    
    double *vel = new double[ ny ];

    for ( int h = 0; h < ( L + 1 ); h++ )
    {
        double d_h = ( double ) h - 0.5;

        vel[ h ] = -4.0 * vx_max / ( L * L ) * ( d_h ) * ( d_h - L );
        
        //cout << "vel[" << h << "] = " << vel[h] << endl;
    }

    for ( int z = 0; z < nz; z++ )
    {
        int h = 0;

        for ( int y = 0; y < ny; y++ )
        {
            int *meio = ini_meio + x + y * nx + z * ny * nx;

            if ( *meio )
            {
                //--------------------------------------------------------------------------------//

                h++;
                
                //cout << "h = " << h << endl;
                
                double vx = vel[ h ];

                //--------------------------------------------------------------------------------//
                
                double *f = ini_f + ( *meio - 1 ) * nvel;

                //------------------- Calcula a soma das incógnitas ------------------------------//

                double sum_f_xminus = f[2] + f[8] + f[10] + f[12] + f[14];

                double sum_f_xzero = f[0] + f[3] + f[4] + f[5] + f[6]
										+ f[15] + f[16] + f[17] + f[18];

                //-------------------- Determinação de vx e imposição de vy e vz -----------------//

                double rho = ( 2.0 * sum_f_xminus + sum_f_xzero ) / ( 1.0 - vx );

                double vy = 0.0;

                double vz = 0.0;

                //--------------------------------------------------------------------------------//

                f[1]  = f[2] + rho * vx / 3.0;

                f[7]  = f[8] + rho * vx / 6.0 + 0.5 * rho * vy - 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[9]  = f[10] + rho * vx / 6.0 - 0.5 * rho * vy + 0.5 * ( f[18] + f[3] + f[16]
                        - f[15] - f[4] - f[17] );

                f[11] = f[12] + rho * vx / 6.0 + 0.5 * rho * vz - 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );

                f[13] = f[14] + rho * vx / 6.0 - 0.5 * rho * vz + 0.5 * ( f[17] + f[5] + f[16]
                        - f[15] - f[6] - f[18] );
            }
        }
    }

    delete[] vel;
}
//================================================================================================//




//===================== Calcula a permeabilidade intrínseca ======================================//
//
//      Input: Quantidade de movimento do fluido (soma_mx), força total sobre o fluido (Q_lost),
//             dimensão do pixel (ftesc), porosidade (phi), viscosidade, dimensão x
//      Output: Permeabilidade
//
//================================================================================================//

double intrinsic_permeability (double soma_mx, double Q_lost, double ftesc, double phi, double visc,
                              int ptos_meio, int nx, int passo, int D_caract )
{

    double mx_med = soma_mx / (double)( nx );

    //------------------------------------------------------------------------------------//

    double Q_lost_med = Q_lost / (double)( nx );
    //double Q_lost_med = sum_force / ( double ) ( nx );

    //------------------------------------------------------------------------------------//

    //phi = 1.0; // Teste
    double k_mts = (ftesc)*(ftesc) * phi * visc * mx_med / Q_lost_med;

    //------------------------------------------------------------------------------------//

    // converte m^2 para mDarcy
    double k_darcy = 10000000.0 * k_mts / ( 0.0000000098697 );

    cout << "\rStep : " << passo << "    k = " << k_darcy << "mDa;    k = " << k_mts
         << "m^2" << endl << endl;

    //------------------------------------------------------------------------------------//

    double Re = ( D_caract ) * ( mx_med /  ptos_meio ) / visc;

    cout << "Reynolds = " << Re << endl << endl;

    ofstream fperme ("k_darcy.dat",ios::app);
    fperme << passo << " " << k_darcy << endl;
    fperme.close();

    ofstream fkm ("k_mts.dat",ios::app);
    fkm << passo << " " << k_mts << endl;
    fkm.close();
    
    return k_darcy;
}

//================================================================================================//




//===================== Lê o arquivo de inicialização (para dois fluidos) ========================//
//
//      Input:
//      Output: nome do arquivo de geometria, tam. do pixel, passos, files, tempo de relaxação,
//              densidade inicial
//
//================================================================================================//

void read_data_twofluid ( char *nome_geo, double *ftesc, int *npassos, int *files, double *tau_R,
                          double *visc_R, double *tau_B, double *visc_B, double *tau_m,
                          double *ro_ini_R, double *ro_ini_B, double *fat_int, double *fat_rec,
                          double *b, double *molhab_R, string* file_rec )

{
     //--------------------------------------------------------------------------------------------//

    char nome_in[15] = "data_in.txt";

    ifstream f_in( nome_in );

    //--------------------------------------------------------------------------------------------//

    char nome_out[15] = "data_out.txt";

    cout << "\nNome do arquivo de dados (saida): " << nome_out << endl;

    //--------------------------------------------------------------------------------------------//

    ofstream fdat( nome_out );

    fdat << "Nome do arquivo de dados (saida): " << nome_out << endl;

    //--------------------------------------------------------------------------------------------//

    f_in >> nome_geo;

    cout << "\nNome do arquivo de geometria: " << nome_geo << endl;

    fdat << "Nome do arquivo de geometria: " << nome_geo << endl;

    //--------------------------------------------------------------------------------------------//

    char st_ft[25];
    f_in >> st_ft;

    f_in >> *ftesc;  // Le a dimensao do dimensão do pixel ( m )

    cout << "\n\nDimensao do pixel = " << *ftesc << endl;

    fdat << "\nDimensao do pixel = " << *ftesc << endl;

    //--------------------------------------------------------------------------------------------//

    char st_steps[25];
    f_in >> st_steps;

    f_in >> *npassos;

    cout << "\nNumero de passos: " << *npassos << endl;

    fdat << "Numero de passos: " << *npassos << endl;

    //--------------------------------------------------------------------------------------------//

    char st_files[15];
    f_in >> st_files;

    f_in >> *files;

    cout << "Numero de arquivos: " << *files << endl;

    //--------------------------------------------------------------------------------------------//

    char st_tau_R[15];
    f_in >> st_tau_R;

    f_in >> *tau_R;

    *visc_R = 1.0 / 3.0 * ( *tau_R - 0.5 );

    cout << "Tempo de relaxacao (Red): " << *tau_R << " => viscosidade R: "
         << *visc_R << endl;

    //--------------------------------------------------------------------------------------------//

    char st_tau_B[15];
    f_in >> st_tau_B;

    f_in >> *tau_B;

    *visc_B = 1.0 / 3.0 * ( *tau_B - 0.5 );

    cout << "Tempo de relaxacao (Blue): " << *tau_B << " => viscosidade B: "
         << *visc_B << endl;

    //--------------------------------------------------------------------------------------------//

    char st_tau_m[15];
    f_in >> st_tau_m;

    f_in >> *tau_m;

    cout << "Tempo de relaxacao (Red_Blue): " << *tau_m << endl;

    //--------------------------------------------------------------------------------------------//

    char st_ro_ini_R[15];
    f_in >> st_ro_ini_R;

    f_in >> *ro_ini_R;

    cout << "Densidade vermelha: " << *ro_ini_R << endl;

    //--------------------------------------------------------------------------------------------//

    char st_ro_ini_B[15];
    f_in >> st_ro_ini_B;

    f_in >> *ro_ini_B;

    cout << "Densidade azul: " << *ro_ini_B << endl;

    //--------------------------------------------------------------------------------------------//

    char st_fat_int[15];
    f_in >> st_fat_int;

    f_in >> *fat_int;

    cout << "Fator de interferencia[0;0.4]: " << *fat_int << endl;

    //--------------------------------------------------------------------------------------------//

    f_in >> st_fat_int;

    f_in >> *fat_rec;

    cout << "Fator de recoloracao[0;1]: " << *fat_rec << endl;

    //--------------------------------------------------------------------------------------------//

    f_in >> st_fat_int;

    f_in >> *b;

    cout << "Termo de exclusao por volume: " << *b << endl;

    //--------------------------------------------------------------------------------------------//

    char st_molhab_R[15];
    f_in >> st_molhab_R;

    f_in >> *molhab_R;

    cout << "Interacao Red_wall [0.0 ; 1.0] : " << *molhab_R << endl;

    //--------------------------------------------------------------------------------------------//
        
    f_in >> *file_rec;
    
    cout << "Inicio da simulação :" << *file_rec << endl;

    f_in.close();

    //--------------------------------------------------------------------------------------------//
}

//================================================================================================//




//===================== Calcula um gradiente de um campo escalar =================================//
//
//      Input: campo escalar, *ini_psi; vetores da rede, endereços do gradiente, posição e tamanho
//      Output: *grad_x, *grad_y, *grad_z
//
//================================================================================================//

void gradient_mirror ( double *ini_psi, double *ini_c, int *ini_meio, double *grad_x, 
					   double *grad_y, double *grad_z, int x, int y, int z, int nx, int ny, int nz, 
					   double *W, double one_over_c_s2 )
{
    double *psi, *c;

    *grad_x = 0.0;
    *grad_y = 0.0;
    *grad_z = 0.0;

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        int  pos_x = x + ( int ) c[0];
        int  pos_y = y + ( int ) c[1];
        int  pos_z = z + ( int ) c[2];
        
        if ( pos_x < 0 ) pos_x = 0;
        if ( pos_x > nx - 1 ) pos_x = nx - 1;
        
        if ( pos_y < 0 ) pos_y = 0;
        if ( pos_y > ny - 1 ) pos_y = ny - 1;
        
        if ( pos_z < 0 ) pos_z = 0;
        if ( pos_z > nz - 1 ) pos_z = nz - 1;

		int *meio = ini_meio + pos_x  +  pos_y * nx + pos_z * nx * ny;
		
		if ( *meio == 0 )
		{
			c = ini_c + i * dim;

			pos_x = ( x - ( int ) c[0] + nx ) % nx;
			pos_y = ( y - ( int ) c[1] + ny ) % ny;
			pos_z = ( z - ( int ) c[2] + nz ) % nz;
		}
		
        psi = ini_psi + pos_x  +  pos_y * nx + pos_z * nx * ny;

        *grad_x = *grad_x + *psi * c[0] * W[i];
        *grad_y = *grad_y + *psi * c[1] * W[i];
        *grad_z = *grad_z + *psi * c[2] * W[i];
    }

}

//================================================================================================//



//===================== Calcula um gradiente de um campo escalar =================================//
//
//      Input: campo escalar, *ini_psi; vetores da rede, endereços do gradiente, posição e tamanho
//      Output: *grad_x, *grad_y, *grad_z
//
//================================================================================================//

void gradient_point ( double *ini_psi, double *ini_c, int *ini_meio, double *grad_x, 
					   double *grad_y, double *grad_z, int x, int y, int z, int nx, int ny, int nz, 
					   double *W, double one_over_c_s2 )
{
    double *psi, *c;

    *grad_x = 0.0;
    *grad_y = 0.0;
    *grad_z = 0.0;

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        int  pos_x = x + ( int ) c[0];
        int  pos_y = y + ( int ) c[1];
        int  pos_z = z + ( int ) c[2];
        
        if ( pos_x < 0 ) pos_x = 0;
        if ( pos_x > nx - 1 ) pos_x = nx - 1;
        
        if ( pos_y < 0 ) pos_y = 0;
        if ( pos_y > ny - 1 ) pos_y = ny - 1;
        
        if ( pos_z < 0 ) pos_z = 0;
        if ( pos_z > nz - 1 ) pos_z = nz - 1;

		int *meio = ini_meio + pos_x  +  pos_y * nx + pos_z * nx * ny;
		
		if ( *meio == 0 )
		{
			c = ini_c + i * dim;

			pos_x = x;
			pos_y = y;
			pos_z = z;
		}
		
        psi = ini_psi + pos_x  +  pos_y * nx + pos_z * nx * ny;

        *grad_x = *grad_x + *psi * c[0] * W[i];
        *grad_y = *grad_y + *psi * c[1] * W[i];
        *grad_z = *grad_z + *psi * c[2] * W[i];
    }

}

//================================================================================================//



//===================== Calcula um gradiente de um campo escalar =================================//
//
//      Input: campo escalar, *ini_psi; vetores da rede, endereços do gradiente, posição e tamanho
//      Output: *grad_x, *grad_y, *grad_z
//
//================================================================================================//

void gradient ( double *ini_psi, double *ini_c, double *grad_x, double *grad_y, double *grad_z,
                int x, int y, int z, int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    double *psi, *c;

    *grad_x = 0.0;
    *grad_y = 0.0;
    *grad_z = 0.0;

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        int  pos_x = ( x + ( int ) c[0] + nx ) % nx;
        int  pos_y = ( y + ( int ) c[1] + ny ) % ny;
        int  pos_z = ( z + ( int ) c[2] + nz ) % nz;

        psi = ini_psi + pos_x  +  pos_y * nx + pos_z * nx * ny;

        *grad_x = *grad_x + *psi * c[0] * W[i];
        *grad_y = *grad_y + *psi * c[1] * W[i];
        *grad_z = *grad_z + *psi * c[2] * W[i];
    }

    *grad_x = *grad_x * one_over_c_s2;
    *grad_y = *grad_y * one_over_c_s2;
    *grad_z = *grad_z * one_over_c_s2;
}

//================================================================================================//



//===================== Impõe a velocidade no ponto de acordo com os vizinhos ====================//
//
//      Input: 
//      Output: 
//
//================================================================================================//

void bound_kr ( int* ini_meio, double *ini_c, double* ini_N, double* ini_N_novo, int x, int y, 
				int z, int nx, int ny, int nz, double *W, double one_over_c_s2 )
{
    double* c;
    
    int* meio;
        
    int pto;
    
    double vx, vy, vz, rho;
    
    int number = 0;

    double sum_vx = 0.0;
    
    double sum_vy = 0.0;
    
    double sum_vz = 0.0;    

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        int  pos_x = ( x + ( int ) c[0] + nx ) % nx;
        int  pos_y = ( y + ( int ) c[1] + ny ) % ny;
        int  pos_z = ( z + ( int ) c[2] + nz ) % nz;

		meio = ini_meio + pos_x + pos_y * nx + pos_z * ny * nx;
		
		if ( *meio )
		{
			pto = ( *meio - 1 );
                        
			double* N = ini_N + ( pto ) * nvel;

			calcula ( N, ini_c, &vx, &vy, &vz, &rho );
			
			sum_vx = sum_vx + vx;
			
			sum_vy = sum_vy + vy;
			
			sum_vz = sum_vz + vz;
			
			number++;
		}        
    }
        
    double vx_avg = sum_vx / ( double ) number;
    
    double vy_avg = sum_vy / ( double ) number;
    
    double vz_avg = sum_vz / ( double ) number; 
    
    meio = ini_meio + x + y * nx + z * ny * nx; 
    
    pto = ( *meio - 1 );
  
    rho = mass ( ini_N + ( pto ) * nvel );    
                        
	double* N_novo = ini_N_novo + ( pto ) * nvel;
        
    dist_eq ( vx_avg, vy_avg, vz_avg, rho, ini_c, N_novo, W, one_over_c_s2 );   
}

//================================================================================================//




//===================== Calcula um gradiente de um campo escalar =================================//
//
//      Input: campo escalar, *ini_psi; vetores da rede, endereços do gradiente, posição e tamanho
//      Output: *grad_x, *grad_y, *grad_z
//
//================================================================================================//

void gradient_with_mirror ( double *ini_psi, double *ini_c, double *grad_x, double *grad_y, 
							double *grad_z, int x, int y, int z, int mirror_x, int mirror_y, int nx, 
							int ny, int nz, double *W, double one_over_c_s2 )
{
    double *psi, *c;

    *grad_x = 0.0;
    *grad_y = 0.0;
    *grad_z = 0.0;

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        int  pos_x = ( x + ( int ) c[0] + nx ) % nx;
        int  pos_y = ( y + ( int ) c[1] + ny ) % ny;
        int  pos_z = ( z + ( int ) c[2] + nz ) % nz;
        
        if ( pos_x > mirror_x )
        {
			pos_x = ( x + ( int ) (- c[0] ) + nx ) % nx;
		}
		
		if ( pos_y > mirror_y )
		{
			pos_y = ( y + ( int ) (- c[1] ) + ny ) % ny;
		}

        psi = ini_psi + pos_x  +  pos_y * nx + pos_z * nx * ny;

        *grad_x = *grad_x + *psi * c[0] * W[i];
        *grad_y = *grad_y + *psi * c[1] * W[i];
        *grad_z = *grad_z + *psi * c[2] * W[i];
    }

    *grad_x = *grad_x * one_over_c_s2;
    *grad_y = *grad_y * one_over_c_s2;
    *grad_z = *grad_z * one_over_c_s2;
}

//================================================================================================//



//===================== Etapa de recoloração (Latva-Koko) ========================================//
//
//      Input: distribution functions, lattice vectors, gradient, magnitude of gradient,
//              mass fractions, density, recolloring factor
//      Output: distribution function
//
//================================================================================================//

void recolloring ( double *N, double *R, double *B, double *ini_c, double mx_m,
                   double my_m, double mz_m, double mod_mM, double concR, double concB, double rho,
                   double beta )

{
    double prod_vm_ci[nvel];

    double cos_phi[nvel];

    double fator = beta * concR * concB;

    double ro_18 = rho / 18.;

    double ro_36 = rho / 36.;

    double *c;

    double c_2;

    //--------------------------------------------------------------------------------------------//

    R[0] = concR * N[0];
    B[0] = concB * N[0];

    for ( int i = 1; i < 7; i++ )
    {
        c = ini_c + i * dim;

        c_2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

        prod_vm_ci[i] = ( mx_m * c[0] + my_m * c[1] + mz_m * c[2] );

        if ( mod_mM ) cos_phi[i] = prod_vm_ci[i] / ( mod_mM * sqrt ( c_2 ) );

        else cos_phi[i] = 0.0;

        R[i] = concR * N[i] - fator * ro_18 * cos_phi[i];
        B[i] = concB * N[i] + fator * ro_18 * cos_phi[i];
    }

    for ( int i = 7; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        c_2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

        prod_vm_ci[i] = ( mx_m * c[0] + my_m * c[1] + mz_m * c[2] );

        if ( mod_mM ) cos_phi[i] = prod_vm_ci[i] / ( mod_mM * sqrt ( c_2 ) );

        else cos_phi[i] = 0.0;

        R[i] = concR * N[i] - fator * ro_36 * cos_phi[i];
        B[i] = concB * N[i] + fator * ro_36 * cos_phi[i];
    }

}
//================================================================================================//




//===================== Imposição de tensão interfacial (Spencer, Halliday & Care) ===============//
//
//      Input: ddistribution function, lattice vectors, factor dependent on the model,
//              unit vectors (normal to the interface)
//      Output: distribution function
//
//================================================================================================//

void imp_interf_tension ( double *N, double *ini_c, double fat, double un_x, double un_y,
                          double un_z, double *W, double one_over_c_s2 )
{
    double *c;

    double c_2;

    double prod_n_ci;

    N[0] = N[0] + fat * W[0] * ( 2./3. );

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        c_2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

        prod_n_ci = ( un_x * c[0] + un_y * c[1] + un_z * c[2] );

        N[i] = N[i] + fat * W[i] * ( prod_n_ci * prod_n_ci - c_2 + 2./3. );
    }
}

//================================================================================================//




//===================== Calcula a "densidade efetiva" ============================================//
//
//      Input: distribution function, densidade inicial
//      Output: densidade efetiva, eff_mass
//
//================================================================================================//

double effetive_mass ( double *f, double rho_0 )
{
    double rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8]
                 + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16]
                 + f[17] + f[18];

    double eff_mass = rho_0 * ( 1.0 - exp ( - rho / rho_0 ) );

    return eff_mass;
}
//================================================================================================//




//===================== Etapa de colisão usando o modelo Shan & Chen =============================//
//
//      Input: distribution function, lattice vectors, velocities, density,
//              acceleration in the x,y,z directions, relaxation time, interaction factor,
//      Output: pos-collisional distribution function
//
//================================================================================================//

void sc_collision ( double *N, double *ini_c, double vx, double vy, double vz, double rho,
                    double acc_x, double acc_y, double acc_z, double tau_N, double G,
                    double *ini_dens_N2, int x, int y, int z, int nx, int ny,
                    int nz )

{
    double cs_2 = 1./3.;

    //double qsi = rho / tau_N;

    double grad_x, grad_y, grad_z;

    gradient ( ini_dens_N2, ini_c, &grad_x, &grad_y, &grad_z, x, y, z, nx, ny, nz, W,
               one_over_c_s2 );

    double vx_eq = vx - tau_N * G * cs_2 * ( grad_x ) + tau_N * acc_x;
    double vy_eq = vy - tau_N * G * cs_2 * ( grad_y ) + tau_N * acc_y;
    double vz_eq = vz - tau_N * G * cs_2 * ( grad_z ) + tau_N * acc_z;

    double feq[nvel];

    dist_eq ( vx_eq, vy_eq, vz_eq, rho, ini_c, feq, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + ( 1.0 / tau_N )*( feq[i] - N[i] );
    }
}

//================================================================================================//




//===================== Etapa de colisão immiscível usando o modelo SFP ==========================//
//
//      Input: distribution functions, lattice vectors, gradient, aceleration,
//                  relaxation times, interaction factor,
//      Output:
//
//================================================================================================//

void sfp_collision ( double *R, double *B, double *ini_c, double mx_m, double my_m,
                     double mz_m, double gxR, double gyR, double gzR, double gxB, double gyB,
                     double gzB, double tau_R, double tau_B, double tau_m, double fat_int,
                     double *W, double one_over_c_s2 )
{

    double *op_col_R = new double[nvel];
    double *op_col_B = new double[nvel];

    double *op_col_RB = new double[nvel];
    double *op_col_BR = new double[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vxR, vyR, vzR, rhoR;
    double vxB, vyB, vzB, rhoB;

    calcula ( R, ini_c, &vxR, &vyR, &vzR, &rhoR );
    calcula ( B, ini_c, &vxB, &vyB, &vzB, &rhoB );

    double rho = rhoR + rhoB;
    double concR = rhoR / rho;
    double concB = 1.0 - concR;

    //------------- Calcula as velocidades modificadas (termo cruzado) ---------------------------//

    double vx_m, vy_m, vz_m;

    vet_unit ( mx_m, my_m, mz_m ,&vx_m, &vy_m, &vz_m );

    double uxR = vxR - fat_int*vx_m;
    double uyR = vyR - fat_int*vy_m;
    double uzR = vzR - fat_int*vz_m;

    double uxB = vxB + fat_int*vx_m;
    double uyB = vyB + fat_int*vy_m;
    double uzB = vzB + fat_int*vz_m;

    //------------ Colisão monofásica e bifásica -------------------------------------------------//

    bgk_op ( R, ini_c, uxB, uyB, uzB, rhoR, gxR, gyR, gzR, tau_m, op_col_RB, W, one_over_c_s2 );
    bgk_op ( B, ini_c, uxR, uyR, uzR, rhoB, gxB, gyB, gzB, tau_m, op_col_BR, W, one_over_c_s2 );

    bgk_op ( R, ini_c, vxR, vyR, vzR, rhoR, gxR, gyR, gzR, tau_R, op_col_R, W, one_over_c_s2 );
    bgk_op ( B, ini_c, vxB, vyB, vzB, rhoB, gxB, gyB, gzB, tau_B, op_col_B, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        R[i] = R[i] + concR * op_col_R[i] + concB * op_col_RB[i];

        B[i] = B[i] + concB * op_col_B[i] + concR * op_col_BR[i];
    }
    
    delete[] op_col_R;
    delete[] op_col_B;
    
    delete[] op_col_RB;
    delete[] op_col_BR;
}

//================================================================================================//




//===================== Colisão immiscível usando o modelo SFP e dois tempos de relaxação ========//
//
//      Input: distribution functions, lattice vectors, gradient, aceleration,
//                  relaxation times, interaction factor,
//      Output:
//
//================================================================================================//

void sfp_collision_TRT ( double *R, double *B, double *ini_c, double mx_m, double my_m,
                     double mz_m, double gxR, double gyR, double gzR, double gxB, double gyB,
                     double gzB, double tau_R_sim, double tau_R_ant, double tau_B_sim, 
                     double tau_B_ant, double tau_m, double fat_int,
                     double *W, double one_over_c_s2 )
{

    double *op_col_R_sim = new double[nvel];
    double *op_col_B_sim = new double[nvel];
	
	double *op_col_R_ant = new double[nvel];
    double *op_col_B_ant = new double[nvel];
    
    double *op_col_RB = new double[nvel];
    double *op_col_BR = new double[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vxR, vyR, vzR, rhoR;
    double vxB, vyB, vzB, rhoB;

    calcula ( R, ini_c, &vxR, &vyR, &vzR, &rhoR );
    calcula ( B, ini_c, &vxB, &vyB, &vzB, &rhoB );

    double rho = rhoR + rhoB;
    double concR = rhoR / rho;
    double concB = 1.0 - concR;

    //------------- Calcula as velocidades modificadas (termo cruzado) ---------------------------//

    double vx_m, vy_m, vz_m;

    vet_unit ( mx_m, my_m, mz_m ,&vx_m, &vy_m, &vz_m );

    double vxR_alt = vxR - fat_int * vx_m;
    double vyR_alt = vyR - fat_int * vy_m;
    double vzR_alt = vzR - fat_int * vz_m;

    double vxB_alt = vxB + fat_int * vx_m;
    double vyB_alt = vyB + fat_int * vy_m;
    double vzB_alt = vzB + fat_int * vz_m;

    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    bgk_op (R, ini_c, vxB_alt,vyB_alt,vzB_alt,rhoR, gxR,gyR,gzR,tau_m, op_col_RB, W, one_over_c_s2);
    bgk_op (B, ini_c, vxR_alt,vyR_alt,vzR_alt,rhoB, gxB,gyB,gzB,tau_m, op_col_BR, W, one_over_c_s2);

    bgk_even ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_sim, op_col_R_sim, W, one_over_c_s2);
    bgk_even ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_sim, op_col_B_sim, W, one_over_c_s2);

	bgk_odd ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_ant, op_col_R_ant, W, one_over_c_s2 );
    bgk_odd ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_ant, op_col_B_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        R[i] = R[i] + concR * ( op_col_R_sim[i] + op_col_R_ant[i] ) + concB * op_col_RB[i];
        B[i] = B[i] + concB * ( op_col_B_sim[i] + op_col_B_ant[i] ) + concR * op_col_BR[i];
    }
    
    delete[] op_col_R_sim;
    delete[] op_col_B_sim;
    
    delete[] op_col_R_ant;
    delete[] op_col_B_ant;
    
    delete[] op_col_RB;
    delete[] op_col_BR;
}

//================================================================================================//



//===================== Colisão immiscível usando o modelo SFP e dois tempos de relaxação ========//
//
//      Input: distribution functions, lattice vectors, gradient, aceleration,
//                  relaxation times, interaction factor,
//      Output:
//
//================================================================================================//

void sfp_collision_TRT_two_factors( double *R, double *B, double *ini_c, double mx_m, double my_m,
                     double mz_m, double gxR, double gyR, double gzR, double gxB, double gyB,
                     double gzB, double tau_R_sim, double tau_R_ant, double tau_B_sim, 
                     double tau_B_ant, double tau_m, double fat_int_R, double fat_int_B,
                     double *W, double one_over_c_s2 )
{

    double *op_col_R_sim = new double[nvel];
    double *op_col_B_sim = new double[nvel];
	
	double *op_col_R_ant = new double[nvel];
    double *op_col_B_ant = new double[nvel];
    
    double *op_col_RB = new double[nvel];
    double *op_col_BR = new double[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vxR, vyR, vzR, rhoR;
    double vxB, vyB, vzB, rhoB;

    calcula ( R, ini_c, &vxR, &vyR, &vzR, &rhoR );
    calcula ( B, ini_c, &vxB, &vyB, &vzB, &rhoB );

    double rho = rhoR + rhoB;
    double concR = rhoR / rho;
    double concB = 1.0 - concR;

    //------------- Calcula as velocidades modificadas (termo cruzado) ---------------------------//

    double vx_m, vy_m, vz_m;

    vet_unit ( mx_m, my_m, mz_m ,&vx_m, &vy_m, &vz_m );

    double vxR_alt = vxR - fat_int_B * vx_m;
    double vyR_alt = vyR - fat_int_B * vy_m;
    double vzR_alt = vzR - fat_int_B * vz_m;

    double vxB_alt = vxB + fat_int_R * vx_m;
    double vyB_alt = vyB + fat_int_R * vy_m;
    double vzB_alt = vzB + fat_int_R * vz_m;

    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    bgk_op (R, ini_c, vxB_alt,vyB_alt,vzB_alt,rhoR, gxR,gyR,gzR,tau_m, op_col_RB, W, one_over_c_s2);
    bgk_op (B, ini_c, vxR_alt,vyR_alt,vzR_alt,rhoB, gxB,gyB,gzB,tau_m, op_col_BR, W, one_over_c_s2);

    bgk_even ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_sim, op_col_R_sim, W, one_over_c_s2);
    bgk_even ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_sim, op_col_B_sim, W, one_over_c_s2);

	bgk_odd ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_ant, op_col_R_ant, W, one_over_c_s2 );
    bgk_odd ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_ant, op_col_B_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        R[i] = R[i] + concR * ( op_col_R_sim[i] + op_col_R_ant[i] ) + concB * op_col_RB[i];
        B[i] = B[i] + concB * ( op_col_B_sim[i] + op_col_B_ant[i] ) + concR * op_col_BR[i];
    }
    
    delete[] op_col_R_sim;
    delete[] op_col_B_sim;
    
    delete[] op_col_R_ant;
    delete[] op_col_B_ant;
    
    delete[] op_col_RB;
    delete[] op_col_BR;
}

//================================================================================================//




//===================== Impõe uma condição inicial aleatória =====================================//
//
//      Input: lattice vectors, geometry, red fluid distributions, blue fluid distributions
//              initial densities, saturation, dimensions
//      Output: red fluid distributions, blue fluid distributions
//
//================================================================================================//

void ini_cond_rand ( double *ini_c, int *ini_meio, double *ini_R, double *ini_B, double rho_ini_R,
                     double rho_ini_B, double sat, int nx, int ny, int nz, double *W,
                     double one_over_c_s2 )

{
    double vx_ini = 0.0;
    double vy_ini = 0.0;
    double vz_ini = 0.0;

    double *R, *B;

    int *meio;

    int saturation = ( int ) ( sat * 1000 );

    //------------------- Inicializa o meio ------------------------------------------------------//

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y * nx + z * nx * ny;

                if ( *meio )
                {
                    R = ini_R + ( *meio - 1 ) * nvel;
                    B = ini_B + ( *meio - 1 ) * nvel;

                    int alea = rand() % 1000;

                    if ( alea < saturation )
                    {
                        dist_eq ( vx_ini, vy_ini, vz_ini, 0.8 * rho_ini_R, ini_c, R, W,
                                  one_over_c_s2 );
                        dist_eq ( vx_ini, vy_ini, vz_ini, 0.2 * rho_ini_B, ini_c, B, W,
                                  one_over_c_s2 );
                    }
                    else
                    {
                        dist_eq ( vx_ini, vy_ini, vz_ini, 0.2 * rho_ini_R, ini_c, R, W,
                                  one_over_c_s2 );
                        dist_eq ( vx_ini, vy_ini, vz_ini, 0.8 * rho_ini_B, ini_c, B, W,
                                  one_over_c_s2 );
                    }
                }
            }
        }
    }
}

//================================================================================================//




//===================== Impõe uma condição inicial via arquivos ==================================//
//
//      Input: lattice vectors, geometry, red fluid distributions, blue fluid distributions
//              initial densities, saturation, initial velocities
//      Output: red fluid distributions, blue fluid distributions
//
//================================================================================================//

void ini_cond_file ( double *ini_c, int *ini_meio, double *ini_R, double *ini_B, double rho_ini_R,
                     double rho_ini_B, string nome_geo, double vx_ini, double vy_ini, double vz_ini,
                     double *W, double one_over_c_s2 )
{   
    int nx, ny, nz;

    ifstream f_ini( nome_geo );

    string line, dump;
	
    stringstream dados;

	for ( int i = 0; i < 4; i++ ) getline( f_ini, dump );
	
	getline( f_ini, line );

    dados << line;    
    
    dados >> dump >> nx >> ny >> nz;
    
    for ( int i = 0; i < 5; i++ ) getline( f_ini, dump );

	dados.clear();

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                int aux;

                f_ini >> aux;

                int *meio;

                meio = ini_meio + x + y*nx + z*ny*nx;

                if ( *meio )
                {
                    double *R, *B;

                    R = ini_R + ( *meio - 1 ) * nvel;

                    B = ini_B + ( *meio - 1 ) * nvel;

                    if ( aux == 2 )
                    {
                        dist_eq ( vx_ini, vy_ini, vz_ini, 0.0,  ini_c, R, W, one_over_c_s2 );

                        dist_eq ( vx_ini, vy_ini, vz_ini, rho_ini_B,  ini_c, B, W, one_over_c_s2 );
                    }
                    if ( aux == 3 )
                    {
                        dist_eq ( vx_ini, vy_ini, vz_ini, rho_ini_R,  ini_c, R, W, one_over_c_s2 );

                        dist_eq ( vx_ini, vy_ini, vz_ini, 0.0,  ini_c, B, W, one_over_c_s2 );
                    }
                }
            }
        }
    }

    f_ini.close();
}

//================================================================================================//





//===================== Retorna o número de pontos de um dos fluidos =============================//
//
//      Input: red fluid distributions, blue fluid distributions
//      Output: number of red fluid sites
//
//================================================================================================//

int n_ptos ( double *ini_F, double *ini_S, int ptos_meio )
{

    double *F;
    double *S;

    int sum = 0;

    for ( int pto = 0; pto < ptos_meio; pto++ )
    {
        F = ini_F + ( pto ) * nvel;
        S = ini_S + ( pto ) * nvel;

        double mass_F = mass ( F );
        double mass_S = mass ( S );

        if ( mass_F > mass_S ) sum++;
    }

    return sum;
}

//================================================================================================//




//===================== Retorna o número de pontos de um dos fluidos (parcial) ===================//
//
//      Input: red fluid distributions, blue fluid distributions, geometry
//      Output: number of red fluid sites
//
//================================================================================================//

int n_ptos_part ( double *ini_F, double *ini_S, int *ini_meio, int x_ini, int x_fim, int y_ini,
                  int y_fim, int z_ini, int z_fim, int nx, int ny, int nz )
{

    double *F;
    double *S;
    int *meio;

    int sum = 0;

    for ( int z = z_ini; z < z_fim; z++ )
    {
        for ( int y = y_ini; y < y_fim; y++ )
        {
            for ( int x = x_ini; x < x_fim; x++ )
            {
                meio = ini_meio + x + y * nx + z * nx * ny;

                if ( *meio )
                {
                    F = ini_F + ( *meio - 1 ) * nvel;
                    S = ini_S + ( *meio - 1 ) * nvel;

                    double mass_F = mass ( F );
                    double mass_S = mass ( S );

                    if ( mass_F > mass_S ) sum++;
                }
            }
        }
    }

    return sum;
}

//================================================================================================//




//===================== Etapa de propagação para a distribuição de mediadores ====================//
//
//      Input: distribution function (*ini_f), new distribution function (*ini_f_new),
//              geometry (*ini_meio), momentum lost (*ini_qlost), wettability, dimensions
//      Output: updated distribution function (*ini_f_new)
//
//================================================================================================//

double propag_med ( int *ini_dir, double *ini_f, double *ini_f_new, int *ini_meio, bool *ini_qlost,
                    double molhab, int ptos_meio )
{

#pragma omp parallel for

   for ( int pto = 0; pto < ptos_meio; pto++ )
    {
        //--------------- Aponta os ponteiros ----------------------------------------------------//
        
		double *f = ini_f + ( pto ) * nvel;

		int *dir = ini_dir + ( pto ) * nvel;

		bool *qlost = ini_qlost + ( pto ) * nvel;
		
		//----------------------------------------------------------------------------------------//

		for ( int i = 0; i < nvel; i++ )
		{
			double *f_new = ini_f_new + dir[i];

			if ( qlost[i] ) *f_new = molhab;

			else *f_new = f[i];
		}
    }
 

    return ( 0 );
}

//================================================================================================//




//===================== Grava os resultados - dois fluidos =======================================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_two_fluid ( int *ini_meio, double *ini_R,  double *ini_B, double *ini_c,
                     unsigned int passo, int nx, int ny, int nz )
{
    double ro;
    double roR;
    double roB;
    double vx, vy, vz;
    double vxR, vyR, vzR;
    double vxB, vyB, vzB;
    double *R;
    double *B;
    int *meio;

    char nomeroR[50];
    char nomeroB[50];
    char nomevel[50];

    sprintf ( nomeroR,"ro_R%06d.vtk", passo );
    sprintf ( nomeroB,"ro_B%06d.vtk", passo );
    sprintf ( nomevel,"vel_%06d.vtk", passo );

    ofstream froR(nomeroR);
    ofstream froB(nomeroB);
    ofstream fvel(nomevel);

    cout << "Gravando.";

    froR << "# vtk DataFile Version 2.0" << endl;
    froR << "Densidade" << endl;
    froR << "ASCII" << endl;
    froR << "DATASET STRUCTURED_POINTS" << endl;
    froR << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    froR << "ASPECT_RATIO 1 1 1" << endl;
    froR << "ORIGIN 0 0 0" << endl;
    froR << "POINT_DATA " << nx*ny*nz << endl;
    froR << "SCALARS densidade double" << endl;
    froR << "LOOKUP_TABLE default" << endl;

    froB << "# vtk DataFile Version 2.0" << endl;
    froB << "Densidade" << endl;
    froB << "ASCII" << endl;
    froB << "DATASET STRUCTURED_POINTS" << endl;
    froB << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    froB << "ASPECT_RATIO 1 1 1" << endl;
    froB << "ORIGIN 0 0 0" << endl;
    froB << "POINT_DATA " << nx*ny*nz << endl;
    froB << "SCALARS densidade double" << endl;
    froB << "LOOKUP_TABLE default" << endl;

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y*nx + z*ny*nx;
                if ( *meio )
                {
                    R = ini_R + ( *meio - 1 ) * nvel;
                    B = ini_B + ( *meio - 1 ) * nvel;

                    calcula ( R, ini_c, &vxR, &vyR, &vzR, &roR );
                    calcula ( B, ini_c, &vxB, &vyB, &vzB, &roB );

                    ro = roR + roB;

                    vx = ( roR * vxR + roB * vxB ) / ro;
                    vy = ( roR * vyR + roB * vyB ) / ro;
                    vz = ( roR * vzR + roB * vzB ) / ro;

                    froR << roR << " ";
                    froB << roB << " ";
                    fvel << vx  << " " << vy << " "<< vz << " ";
                }
                else
                {
                    froR << 0.0 << " ";
                    froB << 0.0 << " ";
                    fvel << 0.0  << " " << 0.0 << " "<< 0.0 << " ";
                }
            }
            froR << endl;
            froB << endl;
            fvel << endl;
        }
    }
    froR.close();
    froB.close();
    fvel.close();

    cout << "\r          " << endl;
}

//================================================================================================//




//===================== Grava os resultados - dois fluidos com exclusão por volume ===============//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_two_fluid_exc ( int *ini_meio, double *ini_R, double *ini_B, double *ini_pv, double *ini_c,
                         unsigned int passo, int nx, int ny, int nz, double *W,
                         double one_over_c_s2 )
{

    double *R;
    double *B;
    int *meio;

    char nomeroR[50];
    char nomeroB[50];
    char nomevel[50];

    sprintf ( nomeroR,"ro_R%06d.vtk", passo );
    sprintf ( nomeroB,"ro_B%06d.vtk", passo );
    sprintf ( nomevel,"vel_%06d.vtk", passo );

    ofstream froR(nomeroR);
    ofstream froB(nomeroB);
    ofstream fvel(nomevel);

    cout << "\nGravando.";

    froR << "# vtk DataFile Version 2.0" << endl;
    froR << "Densidade" << endl;
    froR << "ASCII" << endl;
    froR << "DATASET STRUCTURED_POINTS" << endl;
    froR << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    froR << "ASPECT_RATIO 1 1 1" << endl;
    froR << "ORIGIN 0 0 0" << endl;
    froR << "POINT_DATA " << nx*ny*nz << endl;
    froR << "SCALARS densidade double" << endl;
    froR << "LOOKUP_TABLE default" << endl;

    froB << "# vtk DataFile Version 2.0" << endl;
    froB << "Densidade" << endl;
    froB << "ASCII" << endl;
    froB << "DATASET STRUCTURED_POINTS" << endl;
    froB << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    froB << "ASPECT_RATIO 1 1 1" << endl;
    froB << "ORIGIN 0 0 0" << endl;
    froB << "POINT_DATA " << nx*ny*nz << endl;
    froB << "SCALARS densidade double" << endl;
    froB << "LOOKUP_TABLE default" << endl;

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y*nx + z*ny*nx;

                if ( *meio )
                {
                    R = ini_R + ( *meio - 1 ) * nvel;
                    B = ini_B + ( *meio - 1 ) * nvel;

                    double vxR, vyR, vzR, roR;

                    calcula ( R, ini_c, &vxR, &vyR, &vzR, &roR );

                    double vxB, vyB, vzB, roB;

                    calcula ( B, ini_c, &vxB, &vyB, &vzB, &roB );

                    double rho = roR + roB;

                    double vx_alt = ( roR * vxR + roB * vxB ) / rho;
                    double vy_alt = ( roR * vyR + roB * vyB ) / rho;
                    double vz_alt = ( roR * vzR + roB * vzB ) / rho;

                    //----------------------------------------------------------------------------//

                    double g_x, g_y, g_z;

                    gradient ( ini_pv, ini_c, &g_x, &g_y, &g_z, x, y, z, nx, ny, nz, W,
                               one_over_c_s2 );

                    g_x = -( 1.0 / rho ) * g_x;

                    g_y = -( 1.0 / rho ) * g_y;

                    g_z = -( 1.0 / rho ) * g_z;

                    //----------------------------------------------------------------------------//

                    double vx = vx_alt + 0.5 * g_x;

                    double vy = vy_alt + 0.5 * g_y;

                    double vz = vz_alt + 0.5 * g_z;

                    //----------------------------------------------------------------------------//

                    froR << roR << " ";
                    froB << roB << " ";

                    fvel << vx  << " " << vy << " "<< vz << " ";
                }
                else
                {
                    froR << 0.0 << " ";
                    froB << 0.0 << " ";

                    fvel << 0.0  << " " << 0.0 << " "<< 0.0 << " ";
                }
            }
            froR << endl;
            froB << endl;
            fvel << endl;
        }
    }
    froR.close();
    froB.close();
    fvel.close();

    cout << "\r          " << endl;
}

//================================================================================================//




//===================== Grava os resultados - dois fluidos ( Para o modelo Shan & Chen) ==========//
//
//      Input: geometry, distributions functions(pre and pos), lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_two_fluid_SC ( int *ini_meio, double *ini_R,  double *ini_B, double *ini_R_new,
                        double *ini_B_new, double *ini_c, unsigned int passo, int nx, int ny,
                        int nz )
{
    double *R;
    double *B;
    int *meio;

    char nomeroR[50];
    char nomeroB[50];
    char nomevel[50];

    sprintf ( nomeroR,"ro_R%06d.vtk", passo );
    sprintf ( nomeroB,"ro_B%06d.vtk", passo );
    sprintf ( nomevel,"vel_%06d.vtk", passo );

    ofstream froR(nomeroR);
    ofstream froB(nomeroB);
    ofstream fvel(nomevel);

    cout << "\nGravando.";

    froR << "# vtk DataFile Version 2.0" << endl;
    froR << "Densidade" << endl;
    froR << "ASCII" << endl;
    froR << "DATASET STRUCTURED_POINTS" << endl;
    froR << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    froR << "ASPECT_RATIO 1 1 1" << endl;
    froR << "ORIGIN 0 0 0" << endl;
    froR << "POINT_DATA " << nx*ny*nz << endl;
    froR << "SCALARS densidade double" << endl;
    froR << "LOOKUP_TABLE default" << endl;

    froB << "# vtk DataFile Version 2.0" << endl;
    froB << "Densidade" << endl;
    froB << "ASCII" << endl;
    froB << "DATASET STRUCTURED_POINTS" << endl;
    froB << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    froB << "ASPECT_RATIO 1 1 1" << endl;
    froB << "ORIGIN 0 0 0" << endl;
    froB << "POINT_DATA " << nx*ny*nz << endl;
    froB << "SCALARS densidade double" << endl;
    froB << "LOOKUP_TABLE default" << endl;

    fvel << "# vtk DataFile Version 2.0" << endl;
    fvel << "Velocidade" << endl;
    fvel << "ASCII" << endl;
    fvel << "DATASET STRUCTURED_POINTS" << endl;
    fvel << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    fvel << "ASPECT_RATIO 1 1 1" << endl;
    fvel << "ORIGIN 0 0 0" << endl;
    fvel << "POINT_DATA " << nx*ny*nz << endl;
    fvel << "VECTORS velocidade double" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y*nx + z*ny*nx;
                if ( *meio )
                {
                    R = ini_R + ( *meio - 1 ) * nvel;
                    B = ini_B + ( *meio - 1 ) * nvel;

                    double roR;
                    double roB;

                    double vxR_pre, vyR_pre, vzR_pre;
                    double vxB_pre, vyB_pre, vzB_pre;

                    calcula ( R, ini_c, &vxR_pre, &vyR_pre, &vzR_pre, &roR );
                    calcula ( B, ini_c, &vxB_pre, &vyB_pre, &vzB_pre, &roB );


                    double vxR_pos, vyR_pos, vzR_pos;
                    double vxB_pos, vyB_pos, vzB_pos;

                    calcula ( R, ini_c, &vxR_pos, &vyR_pos, &vzR_pos, &roR );
                    calcula ( B, ini_c, &vxB_pos, &vyB_pos, &vzB_pos, &roB );

                    double rho = roR + roB;

                    double vx, vy, vz;

                    vx = 0.5 * ( roR*vxR_pre + roB*vxB_pre + roR*vxR_pos + roB*vxB_pos ) / rho;
                    vy = 0.5 * ( roR*vyR_pre + roB*vyB_pre + roR*vyR_pos + roB*vyB_pos ) / rho;
                    vz = 0.5 * ( roR*vzR_pre + roB*vzB_pre + roR*vzR_pos + roB*vzB_pos ) / rho;

                    froR << roR << " ";
                    froB << roB << " ";
                    fvel << vx  << " " << vy << " "<< vz << " ";
                }
                else
                {
                    froR << 0.0 << " ";
                    froB << 0.0 << " ";
                    fvel << 0.0  << " " << 0.0 << " "<< 0.0 << " ";
                }
            }
            froR << endl;
            froB << endl;
            fvel << endl;
        }
    }
    froR.close();
    froB.close();
    fvel.close();

    cout << "\r          " << endl;
}

//================================================================================================//





//===================== Grava somente a interface entre dois fluidos =============================//
//
//      Input: geometry, distributions functions, lattice vectors, step, dimensions
//      Output:
//
//================================================================================================//

void rec_interface ( int *ini_meio, double *ini_R, double *ini_B, unsigned int passo, int nx,
                     int ny, int nz )
{

    double *R;
    double *B;
    int *meio;

    string nome = "Interface" + to_string( passo ) + ".vtk";

    ofstream f_inter( nome );

    f_inter << "# vtk DataFile Version 2.0" << endl;
    f_inter << "Interface" << endl;
    f_inter << "ASCII" << endl;
    f_inter << "DATASET STRUCTURED_POINTS" << endl;
    f_inter << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    f_inter << "ASPECT_RATIO 1 1 1" << endl;
    f_inter << "ORIGIN 0 0 0" << endl;
    f_inter << "POINT_DATA " << nx*ny*nz << endl;
    f_inter << "SCALARS densidade double" << endl;
    f_inter << "LOOKUP_TABLE default" << endl;

    for ( int z = 0; z < nz; z++ )
    {
        for ( int y = 0; y < ny; y++ )
        {
            for ( int x = 0; x < nx; x++ )
            {
                meio = ini_meio + x + y*nx + z*ny*nx;

                if ( *meio )
                {
                    R = ini_R + ( *meio - 1 ) * nvel;
                    B = ini_B + ( *meio - 1 ) * nvel;

                    double rhoR = mass( R );

                    double rhoB = mass ( B );

                    double mass_fraction = rhoR / ( rhoR + rhoB );

                    int interface = 1;

                    if ( mass_fraction > 0.4 && mass_fraction < 0.6  ) interface = 0;

                    f_inter << interface << " ";

                }
                else
                {
                    f_inter << 0.0 << " ";
                }
            }
            f_inter << endl;
        }
    }
    f_inter.close();
}

//================================================================================================//




//===================== Etapa de colisão immiscível usando o modelo SHC (Phil. Trans.) ===========//
//
//      Input: distribution functions, lattice vectors, aceleration, relaxation times,
//             interaction factor, recolloring step
//      Output:
//
//================================================================================================//

void shc_coll2011 ( double *R, double *B, double *ini_c, double concR, double concB,
                         double rho, double grad_x, double grad_y, double grad_z, double gx,
                         double gy, double gz, double tau, double fat_int, double beta,
                         double *W, double one_over_c_s2 )
{
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;

    //-------------------- Calcula velocidades e densidades --------------------------------------//

    double prod_conc = concR * concB;

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    //------------------- Colisão monofásica -----------------------------------------------------//

    double *N = new double[nvel];

    for ( int i = 0; i < nvel; i++ ) N[i] = R[i] + B[i];

    bgk_collision ( N, ini_c, tau, gx, gy, gz, W, one_over_c_s2 );

    //------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//

    double fator = ( 1. / tau ) * fat_int * prod_conc * ( 1.5 );

    imp_interf_tension ( N, ini_c, fator, n_x, n_y, n_z, W, one_over_c_s2 );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
    
    delete[] N;
}
//================================================================================================//




//===================== Etapa de colisão immiscível usando o modelo SHC (Phys. Rev.) =============//
//
//      Input: distribution functions, lattice vectors, aceleration, relaxation times,
//             interaction factor, recolloring step
//      Output:
//
//================================================================================================//

void shc_coll2010 ( double *R, double *B, double *ini_c, double concR, double concB,
                         double rho, double grad_x, double grad_y, double grad_z, double gx,
                         double gy, double gz, double tau, double fat_int, double beta,
                         double *W, double one_over_c_s2 )
{
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;

    //-------------------- Calcula velocidades e densidades --------------------------------------//

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    //------------------- Colisão monofásica -----------------------------------------------------//

    double *N = new double[nvel];

    for ( int i = 0; i < nvel; i++ ) N[i] = R[i] + B[i];

    bgk_collision ( N, ini_c, tau, gx, gy, gz, W, one_over_c_s2 );

    //------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//

    double fator = rho * ( 1. / tau ) * fat_int * mod_grad * ( 3.0 );

    imp_interf_tension ( N, ini_c, fator, n_x, n_y, n_z, W, one_over_c_s2 );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
    
    delete[] N;
}
//================================================================================================//




//===================== Colisão immiscível, modelo SHC (Phys. Rev.), Two Rel. Times ==============//
//
//      Input: distribution functions, lattice vectors, aceleration, relaxation times (TRT),
//             interaction factor, recolloring step
//      Output:
//
//================================================================================================//

void shc_coll2010_TRT ( double *R, double *B, double *ini_c, double concR,
                             double concB, double rho, double grad_x, double grad_y, double grad_z,
                             double gx, double gy, double gz, double tau_sim,double tau_ant,
                             double fat_int, double beta, double *W, double one_over_c_s2 )
{
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;

    //-------------------- Calcula velocidades e densidades --------------------------------------//
    
    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    double N[nvel];

    for ( int i = 0; i < nvel; i++ ) N[i] = R[i] + B[i];

    double vx, vy, vz, ro;

    calcula ( N, ini_c, &vx, &vy, &vz, &ro );

    //------------------- Colisão monofásica -----------------------------------------------------//

    double op_col_sim[nvel];

    double op_col_ant[nvel];

    bgk_even ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_sim, op_col_sim, W, one_over_c_s2 );

    bgk_odd ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_ant, op_col_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + op_col_sim[i] + op_col_ant[i];
    }

    //------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//

    double fator =  ( 1. / tau_sim ) * fat_int * mod_grad * ( 3.0 );

    imp_interf_tension ( N, ini_c, fator, n_x, n_y, n_z, W, one_over_c_s2 );

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
}

//================================================================================================//




//===================== Calcula o termo de exclusão por volume para a rede D3Q19 =================//
//
//      Input: gradients, velocities, density, lattice vectors, omega function
//      Output: omega function
//
//================================================================================================//

void omega_exc ( double g_x, double g_y, double g_z, double vx, double vy, double vz, double rho,
                 double *ini_c, double *omega_v, double *W, double one_over_c_s2 )
{

    double *c;

    double v[dim];

    v[0] = vx;
    v[1] = vy;
    v[2] = vz;

    double g[dim];

    g[0] = g_x;
    g[1] = g_y;
    g[2] = g_z;

    double gv = prod_int ( g, v );

    omega_v[0] = - W[0] * rho * one_over_c_s2 *( gv );

    for ( int i = 1; i < nvel; i++ )
    {
        c = ini_c + i * dim;

        double cv = prod_int ( c, v );

        double cg = prod_int ( c, g );

        omega_v[i] = W[i] * rho * one_over_c_s2 *( cg + cg * cv * one_over_c_s2 - gv );
    }
}

//================================================================================================//




//===================== Calcula o termo de exclusão por volume para a rede D3Q19 (model A) =======//
//
//      Input: gradients, velocities, density, lattice vectors, omega function
//      Output: omega function
//
//================================================================================================//

void omega_exc_A ( double pv, double tau, double *ini_c, double *omega_v, double *W, 
					double one_over_c_s2 )
{
    double *c;
    
    double factor = pv * one_over_c_s2 * one_over_c_s2 / ( 2.0 * tau );

    omega_v[0] = factor * W[0] * ( -1.0 );

    for ( int i = 1; i < nvel; i++ ) 
    {
        c = ini_c + i * dim;

        double c2 = prod_int ( c, c );

        omega_v[i] = factor * W[i] * ( c2 - 1.0 );
    }
}

//================================================================================================//




//===================== Colisão, SHC (Phys. Rev.), Two Rel. Times and volume exclusion ===========//
//
//      Input: pressure distribution, distribution functions, lattice vectors, aceleration,
//             relaxation times (TRT), interaction factor, recolloring step
//      Output:
//
//================================================================================================//

void shc_coll2010_TRT_exc ( double *ini_pv, double *R, double *B, double *ini_c,
                            double concR, double concB, double rho, double grad_x, double grad_y,
                            double grad_z, double gx, double gy, double gz, double tau_sim,
                            double tau_ant, double fat_int, double beta,
                            int x, int y, int z, int nx, int ny, int nz, double *W,
                            double one_over_c_s2 )
{
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;

    //-------------------- Calcula velocidades e densidades --------------------------------------//

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    double *N = new double[nvel];

    for ( int i = 0; i < nvel; i++ ) N[i] = R[i] + B[i];

    //--------------------------------------------------------------------------------------------//

    double vx_alt, vy_alt, vz_alt, ro;

    calcula ( N, ini_c, &vx_alt, &vy_alt, &vz_alt, &ro );

    //--------------------------------------------------------------------------------------------//

    double g_x, g_y, g_z;

    gradient ( ini_pv, ini_c, &g_x, &g_y, &g_z, x, y, z, nx, ny, nz, W, one_over_c_s2 );

    g_x = -( 1.0 / rho ) * g_x;

    g_y = -( 1.0 / rho ) * g_y;

    g_z = -( 1.0 / rho ) * g_z;

    //--------------------------------------------------------------------------------------------//

    double vx = vx_alt + 0.5 * g_x;

    double vy = vy_alt + 0.5 * g_y;

    double vz = vz_alt + 0.5 * g_z;

    //------------------- Colisão monofásica -----------------------------------------------------//

    double *op_col_sim = new double[nvel];

    double *op_col_ant = new double[nvel];

    bgk_even ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_sim, op_col_sim, W, one_over_c_s2 );

    bgk_odd ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_ant, op_col_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + op_col_sim[i] + op_col_ant[i];
    }

    //------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//

    double fator = ( 1. / tau_sim ) * fat_int * mod_grad * one_over_c_s2;

    imp_interf_tension ( N, ini_c, fator, n_x, n_y, n_z, W, one_over_c_s2 );

    //------------------- Termo de exclusão por volume -------------------------------------------//

    double *omega_v = new double[nvel];

    omega_exc ( g_x, g_y, g_z, vx, vy, vz, ro, ini_c, omega_v, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + ( 1.0 - 1.0 / ( 2.0 * tau_sim ) ) * omega_v[i];
    }

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
    
    delete[] N;
    
    delete[] op_col_sim;
    delete[] op_col_ant;
    
    delete[] omega_v;
}

//================================================================================================//




//===================== Colisão, SHC (Phys. Rev.), Two Rel. Times and volume exclusion ===========//
//
//      Input: pressure distribution, distribution functions, lattice vectors, aceleration,
//             relaxation times (TRT), interaction factor, recolloring step
//      Output:
//
//================================================================================================//

void shc_coll2010_TRT_exc_E ( double *ini_pv, double *R, double *B, double *ini_c,
                            double concR, double concB, double rho, double grad_x, double grad_y,
                            double grad_z, double gx, double gy, double gz, double tau_sim,
                            double tau_ant, double fat_int, double beta,
                            int x, int y, int z, int nx, int ny, int nz, double *W,
                            double one_over_c_s2 )
{
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;

    //-------------------- Calcula velocidades e densidades --------------------------------------//

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    double *N = new double[nvel];

    for ( int i = 0; i < nvel; i++ ) N[i] = R[i] + B[i];

    //--------------------------------------------------------------------------------------------//

    double vx, vy, vz, ro;

    calcula ( N, ini_c, &vx, &vy, &vz, &ro );

    //--------------------------------------------------------------------------------------------//

    double gp_x, gp_y, gp_z;	// Gradiente da pressão

    gradient ( ini_pv, ini_c, &gp_x, &gp_y, &gp_z, x, y, z, nx, ny, nz, W, one_over_c_s2 );

    gp_x = -( 1.0 / rho ) * gp_x;

    gp_y = -( 1.0 / rho ) * gp_y;

    gp_z = -( 1.0 / rho ) * gp_z;

    //--------------------------------------------------------------------------------------------//
	/*/
    double vx = vx_alt + 0.5 * g_x;

    double vy = vy_alt + 0.5 * g_y;

    double vz = vz_alt + 0.5 * g_z;
	/*/
    //------------------- Colisão monofásica -----------------------------------------------------//

    double *op_col_sim = new double[nvel];

    double *op_col_ant = new double[nvel];

    bgk_even_E ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, gp_x, gp_y, gp_z, tau_sim, op_col_sim, W, 
					one_over_c_s2 );

    bgk_odd_E ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, gp_x, gp_y, gp_z, tau_ant, op_col_ant, W, 
					one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + op_col_sim[i] + op_col_ant[i];
    }

    //------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//

    double fator = ( 1. / tau_sim ) * fat_int * mod_grad * one_over_c_s2;

    imp_interf_tension ( N, ini_c, fator, n_x, n_y, n_z, W, one_over_c_s2 );

    //------------------- Termo de exclusão por volume -------------------------------------------//
	/*/
    double *omega_v = new double[nvel];

    omega_exc ( g_x, g_y, g_z, vx, vy, vz, ro, ini_c, omega_v, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + ( 1.0 - 1.0 / ( 2.0 * tau_sim ) ) * omega_v[i];
    }
	/*/
    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
    
    delete[] N;
    
    delete[] op_col_sim;
    delete[] op_col_ant;
    
    //delete[] omega_v;
}

//================================================================================================//




//===================== Colisão, SHC (Phys. Rev.), Two Rel. Times and volume exclusion ===========//
//
//      Input: pressure distribution, distribution functions, lattice vectors, aceleration,
//             relaxation times (TRT), interaction factor, recolloring step
//      Output:
//
//================================================================================================//

void shc_coll2010_TRT_exc_A ( double *ini_pv, double *R, double *B, double *ini_c,
                            int *ini_meio,
                            double concR, double concB, double rho, double grad_x, double grad_y,
                            double grad_z, double gx, double gy, double gz, double tau_sim,
                            double tau_ant, double fat_int, double beta,
                            int x, int y, int z, int nx, int ny, int nz, double *W,
                            double one_over_c_s2 )
{
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;

    //-------------------- Calcula velocidades e densidades --------------------------------------//

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    double *N = new double[nvel];

    for ( int i = 0; i < nvel; i++ ) N[i] = R[i] + B[i];

    //--------------------------------------------------------------------------------------------//

    double vx, vy, vz, ro;

    calcula ( N, ini_c, &vx, &vy, &vz, &ro );

    //------------------- Colisão monofásica -----------------------------------------------------//

    double *op_col_sim = new double[nvel];

    double *op_col_ant = new double[nvel];

    bgk_even ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_sim, op_col_sim, W, one_over_c_s2 );

    bgk_odd ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_ant, op_col_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + op_col_sim[i] + op_col_ant[i];
    }

    //------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//

    double fator = ( 1. / tau_sim ) * fat_int * mod_grad * one_over_c_s2;

    imp_interf_tension ( N, ini_c, fator, n_x, n_y, n_z, W, one_over_c_s2 );

    //------------------- Termo de exclusão por volume -------------------------------------------//

    double *omega_v = new double[nvel];
    
    double pv = *( ini_pv + x  + y * nx + z * nx * ny );

    omega_exc_A ( pv, tau_sim, ini_c, omega_v, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + omega_v[i];
    }

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
    
    delete[] N;
    
    delete[] op_col_sim;
    delete[] op_col_ant;
    
    delete[] omega_v;
}

//================================================================================================//




//===================== Colisão, SHC (Phys. Rev.), Two Rel. Times and volume exclusion ===========//
//
//      Input: pressure distribution, distribution functions, lattice vectors, aceleration,
//             relaxation times (TRT), interaction factor, recolloring step
//      Output:
//
//================================================================================================//

void shc_coll2010_TRT_exc_mirror ( double *ini_pv, double *R, double *B, double *ini_c,
                            int *ini_meio,
                            double concR, double concB, double rho, double grad_x, double grad_y,
                            double grad_z, double gx, double gy, double gz, double tau_sim,
                            double tau_ant, double fat_int, double beta,
                            int x, int y, int z, int nx, int ny, int nz, double *W,
                            double one_over_c_s2 )
{
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;

    //-------------------- Calcula velocidades e densidades --------------------------------------//

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    double *N = new double[nvel];

    for ( int i = 0; i < nvel; i++ ) N[i] = R[i] + B[i];

    //--------------------------------------------------------------------------------------------//

    double vx_alt, vy_alt, vz_alt, ro;

    calcula ( N, ini_c, &vx_alt, &vy_alt, &vz_alt, &ro );

    //--------------------------------------------------------------------------------------------//

    double g_x, g_y, g_z;

    gradient_mirror ( ini_pv, ini_c, ini_meio, &g_x, &g_y, &g_z, x, y, z, nx, ny, nz, W, 
						one_over_c_s2 );

    g_x = -( 1.0 / rho ) * g_x;

    g_y = -( 1.0 / rho ) * g_y;

    g_z = -( 1.0 / rho ) * g_z;

    //--------------------------------------------------------------------------------------------//

    double vx = vx_alt + 0.5 * g_x;

    double vy = vy_alt + 0.5 * g_y;

    double vz = vz_alt + 0.5 * g_z;

    //------------------- Colisão monofásica -----------------------------------------------------//

    double *op_col_sim = new double[nvel];

    double *op_col_ant = new double[nvel];

    bgk_even ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_sim, op_col_sim, W, one_over_c_s2 );

    bgk_odd ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_ant, op_col_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + op_col_sim[i] + op_col_ant[i];
    }

    //------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//

    double fator = ( 1. / tau_sim ) * fat_int * mod_grad * one_over_c_s2;

    imp_interf_tension ( N, ini_c, fator, n_x, n_y, n_z, W, one_over_c_s2 );

    //------------------- Termo de exclusão por volume -------------------------------------------//

    double *omega_v = new double[nvel];

    omega_exc ( g_x, g_y, g_z, vx, vy, vz, ro, ini_c, omega_v, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + ( 1.0 - 1.0 / ( 2.0 * tau_sim ) ) * omega_v[i];
    }

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
    
    delete[] N;
    
    delete[] op_col_sim;
    delete[] op_col_ant;
    
    delete[] omega_v;
}

//================================================================================================//




//===================== Colisão, SHC (Phys. Rev.), Two Rel. Times and volume exclusion ===========//
//
//      Input: pressure distribution, distribution functions, lattice vectors, aceleration,
//             relaxation times (TRT), interaction factor, recolloring step
//      Output:
//
//================================================================================================//

void shc_coll2010_TRT_exc_point ( double *ini_pv, double *R, double *B, double *ini_c,
                            int *ini_meio,
                            double concR, double concB, double rho, double grad_x, double grad_y,
                            double grad_z, double gx, double gy, double gz, double tau_sim,
                            double tau_ant, double fat_int, double beta,
                            int x, int y, int z, int nx, int ny, int nz, double *W,
                            double one_over_c_s2 )
{
    grad_x = - grad_x;
    grad_y = - grad_y;
    grad_z = - grad_z;

    //-------------------- Calcula velocidades e densidades --------------------------------------//

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    double *N = new double[nvel];

    for ( int i = 0; i < nvel; i++ ) N[i] = R[i] + B[i];

    //--------------------------------------------------------------------------------------------//

    double vx_alt, vy_alt, vz_alt, ro;

    calcula ( N, ini_c, &vx_alt, &vy_alt, &vz_alt, &ro );

    //-------------------- Cálculo de g_alpha para exclusão por volume ---------------------------//

    double gr_pv_x, gr_pv_y, gr_pv_z;

    gradient_point ( ini_pv, ini_c, ini_meio, &gr_pv_x, &gr_pv_y, &gr_pv_z, x, y, z, nx, ny, nz, W, 
						one_over_c_s2 );

    double g_x = -( 1.0 / rho ) * gr_pv_x;

    double g_y = -( 1.0 / rho ) * gr_pv_y;

    double g_z = -( 1.0 / rho ) * gr_pv_z;

    //--------------------------------------------------------------------------------------------//

    double vx = vx_alt + 0.5 * g_x;

    double vy = vy_alt + 0.5 * g_y;

    double vz = vz_alt + 0.5 * g_z;

    //------------------- Colisão monofásica -----------------------------------------------------//

    double *op_col_sim = new double[nvel];

    double *op_col_ant = new double[nvel];

    bgk_even ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_sim, op_col_sim, W, one_over_c_s2 );

    bgk_odd ( N, ini_c, vx, vy, vz, ro, gx, gy, gz, tau_ant, op_col_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + op_col_sim[i] + op_col_ant[i];
    }

    //------------------- Imposição da tensão interfacial (Spencer, Halliday & Care) -------------//

    double fator = ( 1. / tau_sim ) * fat_int * mod_grad * one_over_c_s2;

    imp_interf_tension ( N, ini_c, fator, n_x, n_y, n_z, W, one_over_c_s2 );

    //------------------- Termo de exclusão por volume -------------------------------------------//

    double *omega_v = new double[nvel];

    omega_exc ( g_x, g_y, g_z, vx, vy, vz, ro, ini_c, omega_v, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + ( 1.0 - 1.0 / ( 2.0 * tau_sim ) ) * omega_v[i];
    }

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
    
    delete[] N;
    
    delete[] op_col_sim;
    delete[] op_col_ant;
    
    delete[] omega_v;
}

//================================================================================================//



//===================== Colisão immiscível usando o modelo SFP e dois tempos de relaxação ========//
//
//      Input: distribution functions, lattice vectors, gradient, aceleration,
//                  relaxation times, interaction factor,
//      Output:
//
//================================================================================================//

void sfp_collision_TRT_exc ( double *ini_pv, double *R, double *B, double *ini_c, 
							 double mx_m, double my_m, double mz_m, double gxR, double gyR, 
							 double gzR, double gxB, double gyB, double gzB, double tau_R_sim, 
							 double tau_R_ant, double tau_B_sim, double tau_B_ant, double tau_m, 
							 double fat_int, double *W, double one_over_c_s2, 
							 int x, int y, int z, int nx, int ny, int nz )
{

    double *op_col_R_sim = new double[nvel];
    double *op_col_B_sim = new double[nvel];
	
	double *op_col_R_ant = new double[nvel];
    double *op_col_B_ant = new double[nvel];
    
    double *op_col_RB = new double[nvel];
    double *op_col_BR = new double[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vxR, vyR, vzR, rhoR;
    double vxB, vyB, vzB, rhoB;

    calcula ( R, ini_c, &vxR, &vyR, &vzR, &rhoR );
    calcula ( B, ini_c, &vxB, &vyB, &vzB, &rhoB );

    double rho = rhoR + rhoB;
    double concR = rhoR / rho;
    double concB = 1.0 - concR;
    
    //------------- Altera as velocidades (Exclusão por volume) ----------------------------------//

	double g_x, g_y, g_z;

    gradient ( ini_pv, ini_c, &g_x, &g_y, &g_z, x, y, z, nx, ny, nz, W, one_over_c_s2 );

    g_x = -( 1.0 / rho ) * g_x;

    g_y = -( 1.0 / rho ) * g_y;

    g_z = -( 1.0 / rho ) * g_z;

    //--------------------------------------------------------------------------------------------//

    vxR = vxR + 0.5 * g_x;
    vyR = vyR + 0.5 * g_y;
    vzR = vzR + 0.5 * g_z;

	vxB = vxB + 0.5 * g_x;
    vyB = vyB + 0.5 * g_y;
    vzB = vzB + 0.5 * g_z;

    //------------- Calcula as velocidades modificadas (termo cruzado) ---------------------------//

    double vx_m, vy_m, vz_m;

    vet_unit ( mx_m, my_m, mz_m ,&vx_m, &vy_m, &vz_m );

    double vxR_alt = vxR - fat_int * vx_m;
    double vyR_alt = vyR - fat_int * vy_m;
    double vzR_alt = vzR - fat_int * vz_m;

    double vxB_alt = vxB + fat_int * vx_m;
    double vyB_alt = vyB + fat_int * vy_m;
    double vzB_alt = vzB + fat_int * vz_m;

    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    bgk_op (R, ini_c, vxB_alt,vyB_alt,vzB_alt,rhoR, gxR,gyR,gzR,tau_m, op_col_RB, W, one_over_c_s2);
    bgk_op (B, ini_c, vxR_alt,vyR_alt,vzR_alt,rhoB, gxB,gyB,gzB,tau_m, op_col_BR, W, one_over_c_s2);

    bgk_even ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_sim, op_col_R_sim, W, one_over_c_s2);
    bgk_even ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_sim, op_col_B_sim, W, one_over_c_s2);

	bgk_odd ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_ant, op_col_R_ant, W, one_over_c_s2 );
    bgk_odd ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_ant, op_col_B_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        R[i] = R[i] + concR * ( op_col_R_sim[i] + op_col_R_ant[i] ) + concB * op_col_RB[i];
        B[i] = B[i] + concB * ( op_col_B_sim[i] + op_col_B_ant[i] ) + concR * op_col_BR[i];
    }
    
    //------------------- Termo de exclusão por volume -------------------------------------------//

    double *omega_v_R = new double[nvel];
    double *omega_v_B = new double[nvel];    

    omega_exc ( g_x, g_y, g_z, vxR, vyR, vzR, rhoR, ini_c, omega_v_R, W, one_over_c_s2 );
    omega_exc ( g_x, g_y, g_z, vxB, vyB, vzB, rhoB, ini_c, omega_v_B, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        R[i] = R[i] + ( 1.0 - 1.0 / ( 2.0 * tau_R_sim ) ) * omega_v_R[i];
        B[i] = B[i] + ( 1.0 - 1.0 / ( 2.0 * tau_B_sim ) ) * omega_v_B[i];        
    }
    
    delete[] op_col_R_sim;
    delete[] op_col_B_sim;
    
    delete[] op_col_R_ant;
    delete[] op_col_B_ant;
    
	delete[] op_col_RB;
	delete[] op_col_BR;
	
	delete[] omega_v_R;
	delete[] omega_v_B;
}

//================================================================================================//




//===================== Colisão immiscível usando o modelo SFP e dois tempos de relaxação ========//
//
//      Input: distribution functions, lattice vectors, gradient, aceleration,
//                  relaxation times, interaction factor,
//      Output:
//
//================================================================================================//

void sfp_collision_TRT_exc_point ( double *ini_pv, double *R, double *B, double *ini_c, 
									int *ini_meio, double grad_x, double grad_y, double grad_z, 
									double gxR, double gyR, double gzR, double gxB, double gyB, 
									double gzB, double tau_R_sim, double tau_R_ant, 
									double tau_B_sim, double tau_B_ant, double tau_m, 
									double fat_int, double *W, double one_over_c_s2, 
									int x, int y, int z, int nx, int ny, int nz )
{
	double *op_col_R_sim = new double[nvel];
    double *op_col_B_sim = new double[nvel];
	
	double *op_col_R_ant = new double[nvel];
    double *op_col_B_ant = new double[nvel];
    
    double *op_col_RB = new double[nvel];
    double *op_col_BR = new double[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vxR, vyR, vzR, rhoR;
    double vxB, vyB, vzB, rhoB;

    calcula ( R, ini_c, &vxR, &vyR, &vzR, &rhoR );
    calcula ( B, ini_c, &vxB, &vyB, &vzB, &rhoB );

    double rho = rhoR + rhoB;
    double concR = rhoR / rho;
    double concB = 1.0 - concR;
    
    //------------- Altera as velocidades (Exclusão por volume) ----------------------------------//

	double gr_pv_x, gr_pv_y, gr_pv_z;

    gradient_point ( ini_pv, ini_c, ini_meio, &gr_pv_x, &gr_pv_y, &gr_pv_z, x, y, z, nx, ny, nz, W, 
						one_over_c_s2 );

    double g_x = -( 1.0 / rho ) * gr_pv_x;
    double g_y = -( 1.0 / rho ) * gr_pv_y;
    double g_z = -( 1.0 / rho ) * gr_pv_z;
    
    vxR = vxR + 0.5 * g_x;
    vyR = vyR + 0.5 * g_y;
    vzR = vzR + 0.5 * g_z;
    
    vxB = vxB + 0.5 * g_x;
    vyB = vyB + 0.5 * g_y;
    vzB = vzB + 0.5 * g_z;

    //------------- Calcula as velocidades modificadas (termo cruzado) ---------------------------//

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double vxR_alt = vxR - fat_int * n_x;
    double vyR_alt = vyR - fat_int * n_y;
    double vzR_alt = vzR - fat_int * n_z;

    double vxB_alt = vxB + fat_int * n_x;
    double vyB_alt = vyB + fat_int * n_y;
    double vzB_alt = vzB + fat_int * n_z;

    //------------ Colisão monofásica e bifásica -------------------------------------------------//
        
    bgk_op (R, ini_c, vxB_alt,vyB_alt,vzB_alt,rhoR, gxR,gyR,gzR,tau_m, op_col_RB, W, one_over_c_s2);
    bgk_op (B, ini_c, vxR_alt,vyR_alt,vzR_alt,rhoB, gxB,gyB,gzB,tau_m, op_col_BR, W, one_over_c_s2);

    bgk_even ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_sim, op_col_R_sim, W, one_over_c_s2);
    bgk_even ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_sim, op_col_B_sim, W, one_over_c_s2);

	bgk_odd ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_ant, op_col_R_ant, W, one_over_c_s2 );
    bgk_odd ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_ant, op_col_B_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        R[i] = R[i] + concR * ( op_col_R_sim[i] + op_col_R_ant[i] ) + concB * op_col_RB[i];
        B[i] = B[i] + concB * ( op_col_B_sim[i] + op_col_B_ant[i] ) + concR * op_col_BR[i];
    }
    
    //------------------- Termo de exclusão por volume -------------------------------------------//
    	
	//double vx = concR * vxR + concB * vxB;
	//double vy = concR * vyR + concB * vyB;
	//double vz = concR * vzR + concB * vzB;
		
    double *omega_v_R = new double[nvel];
    double *omega_v_B = new double[nvel];    

    omega_exc ( g_x, g_y, g_z, vxR, vyR, vzR, rhoR, ini_c, omega_v_R, W, one_over_c_s2 );    
    omega_exc ( g_x, g_y, g_z, vxB, vyB, vzB, rhoB, ini_c, omega_v_B, W, one_over_c_s2 );     

    for ( int i = 0; i < nvel; i++ )
    {
        R[i] = R[i] + ( 1.0 - 1.0 * concR / ( 2.0 * tau_R_sim ) ) * omega_v_R[i];       
        B[i] = B[i] + ( 1.0 - 1.0 * concB / ( 2.0 * tau_B_sim ) ) * omega_v_B[i];               
    }
    
    delete[] op_col_R_sim;
    delete[] op_col_B_sim;
    
    delete[] op_col_R_ant;
    delete[] op_col_B_ant;
    
    delete[] op_col_RB;
    delete[] op_col_BR;
    
    delete[] omega_v_R;
    delete[] omega_v_B;
}

//================================================================================================//




//===================== Colisão immiscível usando o modelo SFP e dois tempos de relaxação ========//
//
//      Input: distribution functions, lattice vectors, gradient, aceleration,
//                  relaxation times, interaction factor,
//      Output:
//
//================================================================================================//

void sfp_collision_TRT_exc_point_recoll ( double *ini_pv, double *R, double *B, 
											double *ini_c, int *ini_meio, double grad_x, 
											double grad_y, double grad_z, double gxR, double gyR, 
											double gzR, double gxB, double gyB, 
											double gzB, double tau_R_sim, double tau_R_ant, 
											double tau_B_sim, double tau_B_ant, double tau_m, 
											double fat_int, double beta, double *W, 
											double one_over_c_s2, 
											int x, int y, int z, int nx, int ny, int nz )
{
    double *op_col_R_sim = new double[nvel];
    double *op_col_B_sim = new double[nvel];
	
	double *op_col_R_ant = new double[nvel];
    double *op_col_B_ant = new double[nvel];
    
    double *op_col_RB = new double[nvel];
    double *op_col_BR = new double[nvel];

    //--------------- Calcula vel. e densidade das partículas ------------------------------------//

    double vxR, vyR, vzR, rhoR;
    double vxB, vyB, vzB, rhoB;

    calcula ( R, ini_c, &vxR, &vyR, &vzR, &rhoR );
    calcula ( B, ini_c, &vxB, &vyB, &vzB, &rhoB );

    double rho = rhoR + rhoB;
    double concR = rhoR / rho;
    double concB = 1.0 - concR;
    
    //------------- Altera as velocidades (Exclusão por volume) ----------------------------------//

	double g_x, g_y, g_z;

    gradient_point ( ini_pv, ini_c, ini_meio, &g_x, &g_y, &g_z, x, y, z, nx, ny, nz, W, 
						one_over_c_s2 );

    g_x = -( 1.0 / rho ) * g_x;
    g_y = -( 1.0 / rho ) * g_y;
    g_z = -( 1.0 / rho ) * g_z;

    //--------------------------------------------------------------------------------------------//

    vxR = vxR + 0.5 * g_x;
    vyR = vyR + 0.5 * g_y;
    vzR = vzR + 0.5 * g_z;

	vxB = vxB + 0.5 * g_x;
    vyB = vyB + 0.5 * g_y;
    vzB = vzB + 0.5 * g_z;

    //------------- Calcula as velocidades modificadas (termo cruzado) ---------------------------//

    double n_x, n_y, n_z;

    vet_unit ( grad_x, grad_y, grad_z ,&n_x, &n_y, &n_z );

    double vxR_alt = vxR - fat_int * n_x;
    double vyR_alt = vyR - fat_int * n_y;
    double vzR_alt = vzR - fat_int * n_z;

    double vxB_alt = vxB + fat_int * n_x;
    double vyB_alt = vyB + fat_int * n_y;
    double vzB_alt = vzB + fat_int * n_z;

    //------------ Colisão monofásica e bifásica -------------------------------------------------//
    
    double *N = new double[nvel]; // usado quando há recoloração
    
    bgk_op (R, ini_c, vxB_alt,vyB_alt,vzB_alt,rhoR, gxR,gyR,gzR,tau_m, op_col_RB, W, one_over_c_s2);
    bgk_op (B, ini_c, vxR_alt,vyR_alt,vzR_alt,rhoB, gxB,gyB,gzB,tau_m, op_col_BR, W, one_over_c_s2);

    bgk_even ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_sim, op_col_R_sim, W, one_over_c_s2);
    bgk_even ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_sim, op_col_B_sim, W, one_over_c_s2);

	bgk_odd ( R, ini_c, vxR,vyR,vzR, rhoR, gxR,gyR,gzR, tau_R_ant, op_col_R_ant, W, one_over_c_s2 );
    bgk_odd ( B, ini_c, vxB,vyB,vzB, rhoB, gxB,gyB,gzB, tau_B_ant, op_col_B_ant, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        R[i] = R[i] + concR * ( op_col_R_sim[i] + op_col_R_ant[i] ) + concB * op_col_RB[i];
        B[i] = B[i] + concB * ( op_col_B_sim[i] + op_col_B_ant[i] ) + concR * op_col_BR[i];

		N[i] = R[i] + B[i];
    }
    
    //------------------- Termo de exclusão por volume -------------------------------------------//

    double vx = concR * vxR + concB * vxB;
    double vy = concR * vyR + concB * vyB;
    double vz = concR * vzR + concB * vzB;
    
    double tau_sim = concR * tau_R_sim + concB * tau_B_sim;

    double *omega_v = new double[nvel];

    omega_exc ( g_x, g_y, g_z, vx, vy, vz, rho, ini_c, omega_v, W, one_over_c_s2 );

    for ( int i = 0; i < nvel; i++ )
    {
        N[i] = N[i] + ( 1.0 - 1.0 / ( 2.0 * tau_sim ) ) * omega_v[i];       
    }

    //------------------- Etapa de recoloração (Latva-Koko) --------------------------------------//

	double mod_grad = sqrt ( grad_x * grad_x + grad_y * grad_y + grad_z * grad_z );

    recolloring ( N, R, B, ini_c, grad_x, grad_y, grad_z, mod_grad, concR, concB, rho, beta );
    
    delete[] op_col_R_sim;
    delete[] op_col_B_sim;
    
    delete[] op_col_R_ant;
    delete[] op_col_B_ant;
    
    delete[] op_col_RB;
    delete[] op_col_BR;
}

//================================================================================================//




//===================== Calcula a tensão interfacial pela definição mecânica =====================//
//
//      Input: distribution functions, meio, tamanho x, y, z
//      Output: tensão interfacial
//      Obs: interface sobre o eixo x
//
//================================================================================================//

double calc_inter_tension ( double *ini_R, double *ini_B, int *ini_meio, int nx, int ny, int nz )
{

    int *meio;

    double sum_diff = 0.0;

    int x = nx / 2;

    int z = nz / 2;

    for ( int y = 0; y < ny; y++ )
    {
        meio = ini_meio + x + y * nx + z * ny * nx;

        if ( *meio )
        {
            //----------- Aponta os ponteiros --------------------------------------------//

            double *R, *B;

            R = ini_R + ( *meio - 1 ) * nvel;
            B = ini_B + ( *meio - 1 ) * nvel;

            double Pi_xx = ( R[1] + B[1] + R[2] + B[2] + R[11] + B[11] + R[12] + B[12]
                             + R[13] + B[13] + R[14] + B[14] );

            double Pi_yy = ( R[3] + B[3] + R[4] + B[4] + R[15] + B[15] + R[16] + B[16]
                             + R[17] + B[17] + R[18] + B[18] );
            
          
            double diff = Pi_yy - Pi_xx;

            sum_diff = sum_diff + diff / 2.;
        }
    }
    
    ofstream file_sigma;
    
    file_sigma.open ( "Interf_tension.dat", ios::app );
    
    file_sigma << "Sigma = " << sum_diff << endl;
    
    file_sigma.close();
        
    return sum_diff;

}

//================================================================================================//




//===================== Calcula o campo de velocidades (magnitude) ===============================//
//
//      Input: 
//      Output: 
//
//================================================================================================//

void return_vel_field ( double* ini_R, double* ini_B, double* ini_c, double* ini_vel, int ptos_meio)
{
	double vxR, vyR, vzR, rhoR;
	
	double vxB, vyB, vzB, rhoB;
   
	for ( int pto = 0; pto < ptos_meio; pto++ )
	{
		
        double *R = ini_R + ( pto ) * nvel;
		
		calcula ( R, ini_c, &vxR, &vyR, &vzR, &rhoR );
		
		double *B = ini_B + ( pto ) * nvel;
		
		calcula ( B, ini_c, &vxB, &vyB, &vzB, &rhoB );
		
		double rho = rhoR + rhoB;
		
		double vel_x = ( rhoR * vxR + rhoB * vxB ) / ( rho );
		
		double vel_y = ( rhoR * vyR + rhoB * vyB ) / ( rho );
		
		double vel_z = ( rhoR * vzR + rhoB * vzB ) / ( rho );
		
		double velocity = sqrt( vel_x * vel_x + vel_y * vel_y + vel_z * vel_z );
		
		double *vel = ini_vel + pto;
		
		*vel = velocity;
	}
}

//================================================================================================//




//===================== Compara campos de velocidade =============================================//
//
//      Input: 
//      Output: 
//
//================================================================================================//

double diff_vel_field ( double* ini_vel_old, double* ini_vel, int ptos_meio )
{
   
	double sum_dif_vel2 = 0.0;
	
	double sum_vel2 = 0.0;
   
	for ( int pto = 0; pto < ptos_meio; pto++ )
	{		
        double vel = *( ini_vel + pto );
		
		double vel_0 = *( ini_vel_old + pto );
		
		double dif_vel = vel - vel_0;
		
		double dif_vel2 = dif_vel * dif_vel;
		
		sum_dif_vel2 = sum_dif_vel2 + dif_vel2;
		
		double vel2 = vel * vel;
		
		sum_vel2 = sum_vel2 + vel2;
	}
	
	double delta = sqrt ( sum_dif_vel2 / sum_vel2 );
	
	return delta;
}

//================================================================================================//


