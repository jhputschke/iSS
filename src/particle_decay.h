// Copyright Chun Shen @ 2016
#ifndef SRC_particle_decay_h_
#define SRC_particle_decay_h_

#include <array>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "data_struct.h"
#include "Random.h"

class particle_decay {
 private:
    const std::string table_path_;
    const AfterburnerType afterburner_type_;
    std::vector<particle_info*> resonance_table;
    //! the table entry of a Monte Carlo id (the first one, as a search of
    //! the table from the front would find)
    std::unordered_map<int, particle_info*> resonance_by_pid_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr;

    //! The SMASH list (EOS 91) has 215 channels that are closed at the pole
    //! masses (omega -> pi rho, eta' -> eta sigma, ...): SMASH opens them with
    //! mass-dependent widths. Here a picked channel the pole masses close gets
    //! Breit-Wigner masses for its unstable daughters, and a channel that stays
    //! closed even at their lowest masses is re-drawn among the open ones.
    //! Only for the SMASH list: the UrQMD and PDG lists decay as before.
    const bool smash_list_;
    //! the lightest mass a particle can have and still decay (its own mass if
    //! it is stable): the lower end of its Breit-Wigner as a daughter
    std::unordered_map<int, double> min_mass_cache_;
    //! SMASH list: picked channels closed at the pole masses, channels
    //! re-drawn, mothers left undecayed (no channel open at all)
    long n_closed_at_pole_ = 0;
    long n_redrawn_ = 0;
    long n_undecayed_ = 0;

    particle_info* find_resonance(const int monval) const {
        auto it = resonance_by_pid_.find(monval);
        return (it == resonance_by_pid_.end()) ? nullptr : it->second;
    }

    double min_decay_mass(int monval);
    //! true if the channel is open for a mother of mass M at the daughters'
    //! lowest masses
    bool channel_can_open(double M, const decay_channel_info *channel);
    //! the daughters' masses for a mother of mass M: the pole masses if they
    //! fit, else Breit-Wigner masses for the unstable daughters
    void daughter_masses(double M, const decay_channel_info *channel,
                         std::array<double, 5> &masses);

 public:
    particle_decay(std::shared_ptr<RandomUtil::Random> ran_gen,
                   AfterburnerType afterburner_type,
                   std::string table_path="iSS_tables");
    ~particle_decay();

    //! This function reads in resonance decay table
    int read_resonances_list();

    //! This is a test function to check whether the resonance table is
    //! read in correctly
    void check_resonance_table();

    //! This function returns particle width in GeV
    double get_particle_width(iSS_Hadron *part);

    //! This function checks whether the particle is stable
    int check_particle_stable(iSS_Hadron *part);

    //! This function returns the electric charge of particle
    int get_particle_charge(int monval);

    //! This function returns the baryon number of particle
    int get_particle_baryon_number(int monval);

    //! This function returns the strange number of particle
    int get_particle_strange_number(int monval);

    //! This function returns the particle mass for a given particle id
    double get_particle_mass(int POI_monval);

    //! This is a shell function to perform resonance decays. Returns false
    //! if the mother could not decay (SMASH list only, no channel open at its
    //! mass): then nothing is added, and the caller keeps the mother.
    bool perform_decays(iSS_Hadron *mother,
                        std::vector<iSS_Hadron>* daughter_list);

    //! SMASH list counters (see smash_list_)
    long get_n_closed_at_pole() const { return n_closed_at_pole_; }
    long get_n_redrawn() const { return n_redrawn_; }
    long get_n_undecayed() const { return n_undecayed_; }


    //! This function perform two body decay
    void perform_two_body_decay(iSS_Hadron *mother,
                                iSS_Hadron *daughter1,
                                iSS_Hadron *daughter2);

    //! This function perform 3 body decays
    void perform_three_body_decay(iSS_Hadron *mother,
                                  iSS_Hadron *daughter1,
                                  iSS_Hadron *daughter2,
                                  iSS_Hadron *daughter3);

    //! This function sample mother particle mass according to breit-wigner
    //! distribution
    double sample_breit_wigner(double mass, double width, double M_min);

    //! the same, truncated to M_min <= m <= M_max
    double sample_breit_wigner(double mass, double width, double M_min,
                               double M_max);
};

#endif  // SRC_particle_decay_h_
