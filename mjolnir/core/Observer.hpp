#ifndef MJOLNIR_CORE_OBSERVER
#define MJOLNIR_CORE_OBSERVER
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/progress_bar.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace mjolnir
{


// TODO: consider adding ObserverBase and implement Observers for each format
template<typename traitsT>
class Observer
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef progress_bar<50> progress_bar_type;

  public:

    Observer(const std::string& filename_prefix, bool output_progress = false)
      : output_progress_(output_progress), progress_bar_(1), prefix_(filename_prefix),
        xyz_name_(filename_prefix + std::string(".xyz")),
        vel_name_(filename_prefix + std::string(".vel")),
        ene_name_(filename_prefix + std::string(".ene"))
    {
        // clear the contents
        {
            std::ofstream ofs(this->xyz_name_);
            if(not ofs.good())
            {
                throw std::runtime_error("file open error: " + this->xyz_name_);
            }
            ofs.close();
        }
        {
            std::ofstream ofs(this->vel_name_);
            if(not ofs.good())
            {
                throw std::runtime_error("file open error: " + this->vel_name_);
            }
            ofs.close();
        }
        {
            std::ofstream ofs(this->ene_name_);
            if(not ofs.good())
            {
                throw std::runtime_error("file open error: " + this->ene_name_);
            }
            ofs.close();
        }
    }
    ~Observer() = default;

    void initialize(const system_type& sys, const forcefield_type& ff,
                    const std::size_t total_step)
    {
        this->progress_bar_.reset(total_step); // set total_step

        std::ofstream ofs(this->ene_name_, std::ios::app);
        ofs << "# timestep  " << ff.list_energy_name() << " kinetic_energy\n";
        return;
    }

    void output(const std::size_t step,
                const system_type& sys, const forcefield_type& ff);

    void output_progress   (const std::size_t step);
    void output_progress_LF(const std::size_t step);

    std::string const& prefix() const noexcept {return prefix_;}

  private:

    real_type calc_kinetic_energy(const system_type& sys) const
    {
        real_type k = 0.0;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            k += length_sq(sys[i].velocity) * sys[i].mass;
        }
        return k * 0.5;
    }

  private:

    bool output_progress_;
    std::string prefix_;
    std::string xyz_name_;
    std::string vel_name_;
    std::string ene_name_;
    progress_bar_type progress_bar_;
};

template<typename traitsT>
inline void Observer<traitsT>::output_progress(const std::size_t step)
{
    // XXX consider introducing template argument to remove this if-branching
    //     at the compile time
    if(this->output_progress_)
    {
        std::cerr << progress_bar_.format(step);
    }
    return;
}
template<typename traitsT>
inline void Observer<traitsT>::output_progress_LF(const std::size_t step)
{
    // XXX consider introducing template argument to remove this if-branching
    //     at the compile time
    if(this->output_progress_)
    {
        std::cerr << progress_bar_.format(step) << std::endl;
    }
    return;
}

template<typename traitsT>
inline void Observer<traitsT>::output(
    const std::size_t step, const system_type& sys, const forcefield_type& ff)
{
    std::ofstream ofs(xyz_name_, std::ios::app);
    ofs << sys.size() << "\nstep = " << step << '\n';
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& p = sys.position(i);
        ofs << sys.name(i) << ' ' << std::fixed << std::setprecision(8)
            << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
    }
    ofs.close();

    ofs.open(vel_name_, std::ios::app);
    ofs << sys.size() << "\nstep = " << step << '\n';
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& v = sys.velocity(i);
        ofs << sys.name(i) << ' ' << std::fixed << std::setprecision(8)
            << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
    }
    ofs.close();

    ofs.open(ene_name_, std::ios::app);
    // if the width exceeds, operator<<(std::ostream, std::string) ignores
    // ostream::width and outputs whole string.
    ofs << std::setw(11) << std::left << std::to_string(step) << ' ';
    ofs << ff.dump_energy(sys) << ' ';
    ofs << std::setw(14) << std::right << this->calc_kinetic_energy(sys) << '\n';
    ofs.close();

    return ;
}

} // mjolnir
#endif /* MJOLNIR_CORE_OBSERVER */
