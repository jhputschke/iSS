// Copyright @ 2020 Chun Shen
// The Fastest Spectra Sampler in the West

#ifndef SRC_FSSW_H_
#define SRC_FSSW_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <memory>

#include "Table.h"
#include "TableFunction.h"
#include "ParameterReader.h"
#include "particle_decay.h"
#include "Spectators.h"
#include "data_struct.h"
#include "Random.h"
#include "pretty_ostream.h"
#include "MomentumSamplerShell.h"
#include "eos_4D_foriSS.h"

class FSSW {
 private:
    int hydro_mode;   // switch for (2+1)-d or (3+1)-d hypersurface
    int flag_PCE_;
    int echoLevel_;

    pretty_ostream messager_;

    const std::string path_;
    const std::string table_path_;
    const AfterburnerType afterburner_type_;

    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;
    std::shared_ptr<MomentumSamplerShell> momentum_sampler_ptr_;

    int number_of_repeated_sampling;
    int maximumSamplingEvents_;
    int local_charge_conservation;

    ParameterReader *paraRdr; // used to pass-in parameters
    bool bUSE_OSCAR_FORMAT;
    bool bUSE_OSCAR2013;
    bool bUSE_GZIP_FORMAT;
    bool bUSE_BINARY_FORMAT;
    bool bINCLUDE_SHEAR_DELTAF, bINCLUDE_BULK_DELTAF, bINCLUDE_DIFFUSION_DELTAF;
    int NEoS_deltaf_kind_;        // 0: 14-mom, 1: CE
    int bulk_deltaf_kind_;
    int iEoS_MUSIC_;
    bool flagRegVisYield_;

    // dN/(dxt dy) for one particle species
    std::vector<double> dN_dxtdy_for_one_particle_species;

    // delta f coefficients of every cell, computed once per surface: they
    // depend on the cell only, not on the species or the sampled hadron.
    // cellVisCoeffs_[l*kMaxVisCoeffs_ + i], cellVisCoeffsSize_[l] of them.
    static const int kMaxVisCoeffs_ = 8;
    bool cellVisCoeffsReady_ = false;
    std::vector<double> cellVisCoeffs_;
    std::vector<unsigned char> cellVisCoeffsSize_;

    // what the yield loop reads of every cell, packed once per surface
    // (prepare_cell_yield_inputs): the loop streams every cell once per
    // block of species and is bound by memory traffic, and FO_surf_LRF is
    // 168 bytes per cell for the 24 it needs
    bool cellYieldInputsReady_ = false;
    std::vector<float> cell_dsigma_u_, cell_T_;
    std::vector<float> cell_muB_, cell_muS_, cell_muQ_, cell_bulkPi_;
    // diffusion delta f only
    std::vector<double> cell_dsigma_q_, cell_qmu_coeff_, cell_prefactor_qmu_;

    // the yields of a block of species from one pass over the cells,
    // dN_block_[s*FO_length + l]
    static const int kSpeciesBlock_ = 8;
    std::vector<double> dN_block_;

    // the inverse CDF over the cells, reused from species to species
    std::vector<double> invCDF_storage_;

    // (unused: the decay chains are per thread now, SamplerWorker; kept so
    // that the members below keep their offsets for code built against the
    // earlier header)
    std::vector<iSS_Hadron> decay_chain_, decay_daughters_;

    // correlated sampling (correlated_sampling = 1): the random numbers are
    // addressed by the physical cell, not by the order of the draws, so two
    // surfaces that agree in a region give the same hadrons there. Cells are
    // grouped into blocks of quantized (tau, x, y, eta); within a block they
    // keep the surface's order.
    bool correlated_sampling_ = false;
    double corr_block_dtau_ = 0.5, corr_block_dx_ = 1.0, corr_block_deta_ = 0.5;
    std::vector<long> block_order_;   // cell indices sorted by (block, index)
    std::vector<long> block_start_;   // block b: block_order_[start[b], start[b+1])
    std::vector<std::array<uint32_t, 2>> block_word_;  // block b's address

    int number_of_chosen_particles;

    // store particle index;
    // the sampling process follows the order specified by this table
    std::vector<int> chosen_particles_sampling_table;

    // list for information for all particles
    std::vector<particle_info> particles;

    // list for information for all fluid cells
    long FO_length;
    const std::vector<FO_surf_LRF> &FOsurf_ptr;

    inline long determine_number_to_sample(
                    double dN, int model=1, double para1=0);

    const gsl_rng_type *gsl_type_random_number;
    gsl_rng *gsl_random_r;

    bool particles_are_the_same(int, int);

    //array for bulk delta f coefficients
    std::unique_ptr<Table> bulkdf_coeff_;

    // table parameter for diffusion deltaf coefficient
    int deltaf_qmu_coeff_table_length_T;
    int deltaf_qmu_coeff_table_length_mu;
    double delta_qmu_coeff_table_T0, delta_qmu_coeff_table_mu0;
    double delta_qmu_coeff_table_dT, delta_qmu_coeff_table_dmu;
    double **deltaf_qmu_coeff_tb;       // kappa_B [1/fm^3]

    // table parameter for bulk deltaf coefficient 14 moments
    int deltaf_bulk_coeff_14mom_table_length_T_;
    int deltaf_bulk_coeff_14mom_table_length_mu_;
    double deltaf_bulk_coeff_14mom_table_T0_;
    double deltaf_bulk_coeff_14mom_table_mu0_;
    double deltaf_bulk_coeff_14mom_table_dT_;
    double deltaf_bulk_coeff_14mom_table_dmu_;
    double **deltaf_bulk_coeff_14mom_c0_tb_;
    double **deltaf_bulk_coeff_14mom_c1_tb_;
    double **deltaf_bulk_coeff_14mom_c2_tb_;

    EOS_4D_foriSS eos_4d_;

    // table parameter for bulk deltaf coefficient CE for NEOS BQS
    int deltaf_coeff_NEOSBQS_table_length_e_;
    int deltaf_coeff_NEOSBQS_table_length_nB_;
    double deltaf_coeff_NEOSBQS_table_de_;
    std::vector<std::vector<double>> deltaf_coeff_22mom_NEOSBQS_tb_;
    std::vector<std::vector<double>> deltaf_coeff_CE_NEOSBQS_tb_;

    // arrays to speed up computing particle yield
    double sf_dx, sf_x_min, sf_x_max;
    int sf_tb_length;
    int sf_expint_truncate_order;
    double** sf_bessel_Kn;
    double** sf_expint_En;

    int flag_output_samples_into_files;
    int flag_store_samples_in_memory;

    //! particle decay
    bool flag_perform_decays_;
    std::unique_ptr<particle_decay> decayer_ptr_;
    //! decay_smash_list = 1: iSS decays the SMASH hadron list (EOS 91) too,
    //! for runs without SMASH. 0 (default): the SMASH list is left to SMASH,
    //! as before.
    bool decay_smash_list_ = false;
    //! decays happen: perform_decays, and not a SMASH list left to SMASH
    bool decays_on() const {
        return flag_perform_decays_
               && (afterburner_type_ != AfterburnerType::SMASH
                   || decay_smash_list_);
    }

    //! keep_primaries = N: the undecayed hadrons of the first N samples, as
    //! sampled (before any decay), next to Hadron_list: for an afterburner
    //! run later. Takes no random numbers, so Hadron_list is the same.
    int keep_primaries_ = 0;
    std::vector<std::vector<iSS_Hadron>> primary_list_;

    //! spectators
    bool flag_spectators_;
    std::unique_ptr<Spectators> spectators_ptr_;

    std::vector< std::vector<iSS_Hadron>* >* Hadron_list;

    std::string OSCAR_header_filename_, OSCAR_output_filename_;

 public:
    FSSW(std::shared_ptr<RandomUtil::Random> ran_gen,
           Table* chosen_particles_in,
           std::vector<particle_info> particles_in,
           const std::vector<FO_surf_LRF> &FOsurf_ptr_in, int flag_PCE,
           ParameterReader* paraRdr_in, string path, string table_path,
           AfterburnerType afterburner_type);
    ~FSSW();


    void initialize_special_function_arrays();
    double get_special_function_K1(double arg);
    double get_special_function_K2(double arg);
    double get_special_function_K3(double arg);
    void get_special_function_En(double arg, std::vector<double> &results);

    void shell();  // it all starts here...

    void combine_samples_to_OSCAR();
    void combine_samples_to_OSCAR2013();
    void combine_samples_to_gzip_file();
    void combine_samples_to_binary_file();

    void calculate_dN_analytic(const particle_info* particle, double mu,
                               double Temperature,
                               std::array<double, 6> &results);

    // the following variables need to be set first in order to
    // call this function
    void calculate_dN_dxtdy_for_one_particle_species(const int particle_idx);
    void calculate_dN_dxtdy_for_species_block(
                const std::vector<int> &species_idx,
                double *__restrict__ dN_out,
                double *block_yield_out = nullptr);
    void sample_using_dN_dxtdy_4all_particles_conventional();
    void sample_one_species_conventional(const int real_particle_idx,
                                         const double *dN_cell);

    void prepare_cell_blocks();
    void release_cell_blocks();
    void sample_one_species_correlated(const int real_particle_idx,
                                       const double *dN_cell,
                                       const double *block_yields);

    //! Correlated sampling on several threads (OpenMP). Every draw is
    //! addressed by (seed, species, block, hadron) in a Philox stream, so a
    //! block's hadrons do not depend on the thread that samples it. The
    //! (species, chunk of blocks) items of a group of species are shared out
    //! among the threads; each thread appends to a buffer of its own, and
    //! the items are then appended to the samples in (species, block) order,
    //! as one thread would have: the output is the same for any number of
    //! threads. Each thread has its own generator, momentum sampler and
    //! decayer (they hold the generator and caches) and scratch vectors;
    //! worker 0 uses the sampler's own.

    //! where a block's hadrons go: straight into the samples (one thread),
    //! or into a thread's buffer (several threads)
    struct BlockSink {
        std::vector<std::vector<iSS_Hadron>*> *samples = nullptr;
        std::vector<std::vector<iSS_Hadron>> *primaries = nullptr;
        std::vector<std::pair<int, iSS_Hadron>> hadrons, kept;
        void hadron(const int isample, const iSS_Hadron &h) {
            if (samples) (*samples)[isample]->push_back(h);
            else hadrons.emplace_back(isample, h);
        }
        void primary(const int isample, const iSS_Hadron &h) {
            if (primaries) (*primaries)[isample].push_back(h);
            else kept.emplace_back(isample, h);
        }
    };
    struct SamplerWorker {
        std::shared_ptr<RandomUtil::Random> ran_gen;
        std::shared_ptr<MomentumSamplerShell> momentum_sampler;
        particle_decay *decayer = nullptr;
        std::unique_ptr<particle_decay> own_decayer;
        std::vector<iSS_Hadron> decay_chain, decay_daughters;
        std::vector<double> visCoeffs_local, block_cdf;
        BlockSink sink;   // this thread's buffer
    };
    //! one (species, chunk of blocks) item: whose buffer, which entries
    struct ItemOutput {
        int thread;
        size_t h0, h1, k0, k1;
    };
    static const long kBlocksPerItem_ = 64;
    void prepare_workers();
    void release_workers();
    void sample_species_group_correlated(const std::vector<int> &species_idx,
                                         const double *dN_group);
    void sample_block_correlated(SamplerWorker &w, BlockSink &sink,
                                 const particle_info *particle,
                                 const double *dN_cell,
                                 const double block_yield, const long b,
                                 const uint32_t seed, const int n_keep,
                                 const bool decays, const double y_LB,
                                 const double y_RB, const double y_factor);
    void decay_one_primary(SamplerWorker &w, const iSS_Hadron &primary,
                           const int isample, BlockSink &sink);

    void getCellVisCoefficients(const FO_surf_LRF *surf,
                                std::vector<double> &visCoefficients);
    void prepare_cell_visCoefficients();
    void release_cell_visCoefficients();
    void prepare_cell_yield_inputs();
    void release_cell_yield_inputs();

    void getbulkvisCoefficients(const double Tdec,
                                std::vector<double> &bulkvisCoefficients);
    void getbulkvisCoefficients(const double Tdec, const double mu_B,
                                std::vector<double> &bulkvisCoefficients);
    void getCENEOSBQSCoefficients(const double Edec, const double nB,
                                  std::vector<double> &visCoefficients);
    void get22momNEOSBQSCoefficients(const double Edec, const double nB,
                                     std::vector<double> &visCoefficients);

    void load_bulk_deltaf_14mom_table(string filepath);
    void load_CE_deltaf_NEOSBQS_table(string filepath);
    void load_22mom_deltaf_NEOSBQS_table(string filepath);
    void load_deltaf_qmu_coeff_table(std::string filename);
    double get_deltaf_qmu_coeff(double T, double muB);

    void check_samples_in_memory();
    int get_number_of_sampled_events() {return(Hadron_list->size());};
    int get_number_of_particles(int iev) {
        return((*Hadron_list)[iev]->size());
    };

    iSS_Hadron get_hadron(int iev, int ipart) {
        return((*(*Hadron_list)[iev])[ipart]);
    };

    std::vector<iSS_Hadron>* get_hadron_list_iev(const int iev) {
        return(Hadron_list->at(iev));
    }

    //! keep_primaries: the samples with primaries (the first N), and one
    int get_number_of_primary_samples() const {
        return(static_cast<int>(primary_list_.size()));
    }
    std::vector<iSS_Hadron>* get_primary_list_iev(const int iev) {
        return(&primary_list_.at(iev));
    }

    void addSpectatorsToHadronList();

    void computeAvgTotalEnergyMomentum();

    void perform_resonance_feed_down(
                std::vector< std::vector<iSS_Hadron>* >* input_particle_list);
    int compute_number_of_sampling_needed(int number_of_particles_needed);
    double get_deltaf_bulk(
        const double mass, const double pdotu, const double bulkPi,
        const double Tdec, const int sign, const int baryon,
        const int strange, const int charge,
        const double f0, const double *bulkvisCoefficients);
    int sample_momemtum_from_a_fluid_cell(
        const double mass, const int sign,
        const int baryon, const int strange, const int charge,
        const FO_surf_LRF *surf,
        const double *visCoefficients,
        const double deltaf_qmu_coeff,
        double &pT, double &phi, double &y_minus_eta_s) {
        return sample_momemtum_from_a_fluid_cell(
            *ran_gen_ptr, *momentum_sampler_ptr_, mass, sign, baryon,
            strange, charge, surf, visCoefficients, deltaf_qmu_coeff,
            pT, phi, y_minus_eta_s);
    }
    //! the same with a given generator and momentum sampler (one thread's)
    int sample_momemtum_from_a_fluid_cell(
        RandomUtil::Random &ran_gen,
        const MomentumSamplerShell &momentum_sampler,
        const double mass, const int sign,
        const int baryon, const int strange, const int charge,
        const FO_surf_LRF *surf,
        const double *visCoefficients,
        const double deltaf_qmu_coeff,
        double &pT, double &phi, double &y_minus_eta_s);
    iSS_Hadron make_sampled_particle(
        const FO_surf_LRF *surf, const int particle_monval, const double mass,
        const double pT, const double phi,
        const double y_minus_eta_s, const double eta_s) const;
    void add_one_sampled_particle(
        const int repeated_sampling_idx, const FO_surf_LRF *surf,
        const int particle_monval, const double mass,
        const double pT, const double phi,
        const double y_minus_eta_s, const double eta_s);
    void boost_vector_back_to_lab_frame(iSS_data::VecD4 &p_LRF,
                                        iSS_data::VecD4 &p_lab,
                                        iSS_data::VecD4 &umu) const;
    double bilinearInterp(std::vector<std::vector<double>>&mat,
                          int idx_e1, int idx_e2, int idx_nB1, int idx_nB2,
                          double x_fraction, double y_fraction);
    int getParticleInfo(const int monval) const;

 private:
    // (last, so that the members above keep their offsets for code built
    // against the earlier header)
    std::vector<std::unique_ptr<SamplerWorker>> workers_;
    std::vector<ItemOutput> item_output_;
    //! correlated sampling: the yield of every (species s of the current
    //! group, block b) at [s*n_blocks + b], summed in the cells' order
    std::vector<double> block_yield_;
    //! correlated sampling: the per-cell arrays (delta f coefficients, yield
    //! inputs, yields) are kept in block order, position j holding cell
    //! cell_order_[j] (block_order_), so that a block's cells are
    //! contiguous; nullptr: in the surface's order
    const long *cell_order_ = nullptr;
    long cell_at(const long l) const {
        return cell_order_ ? cell_order_[l] : l;
    }
};

#endif  // SRC_FSSW_H_
