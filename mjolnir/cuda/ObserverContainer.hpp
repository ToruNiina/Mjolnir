#ifndef MJOLNIR_CUDA_OBSERVER_CONTAINER_HPP
#define MJOLNIR_CUDA_OBSERVER_CONTAINER_HPP
#include <mjolnir/core/ObserverContainer.hpp>
#include <mjolnir/cuda/CUDASimulatorTraits.hpp>
#include <mjolnir/cuda/System.hpp>

// This class manages several different XXXObservers to output
// positions, energies, topologies, and others.
//
// As mentioned at the top of cuda/System.hpp, this CUDA implementation does
// all the calculations on device. That means that only observers needs to
// pull system configuration back to CPU.
namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
class ObserverContainer<CUDASimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type        = CUDASimulatorTraits<realT, boundaryT>;
    using observer_base_type = ObserverBase<traits_type>;
    using observer_base_ptr  = std::unique_ptr<observer_base_type>;
    using real_type          = typename observer_base_type::real_type;
    using coordinate_type    = typename observer_base_type::coordinate_type;
    using system_type        = typename observer_base_type::system_type;
    using forcefield_type    = typename observer_base_type::forcefield_type;
    using progress_bar_type  = progress_bar</*width of bar = */50>;

  public:

    explicit ObserverContainer(bool output_progress = false)
        : output_progress_(output_progress && io::detail::isatty(std::cerr))
    {}
    ~ObserverContainer() = default;

    ObserverContainer(const ObserverContainer&) = default;
    ObserverContainer(ObserverContainer&&)      = default;
    ObserverContainer& operator=(const ObserverContainer&) = default;
    ObserverContainer& operator=(ObserverContainer&&)      = default;

    // open files, write header and so on.
    void initialize(const std::size_t total_step, const real_type dt,
                    const system_type& sys, const forcefield_type& ff)
    {
        for(const auto& obs : observers_)
        {
            obs->initialize(total_step, dt, sys, ff);
        }

        this->progress_bar_.reset(total_step); // set total_step as 100%.
    }
    // call if system or forcefield is changed.
    void update(const std::size_t step, const real_type dt,
                const system_type& sys, const forcefield_type& ff)
    {
        for(const auto& obs : this->observers_)
        {
            obs->update(step, dt, sys, ff);
        }
    }
    // output the current state.
    void output(const std::size_t step, const real_type dt,
                system_type&      sys_, const forcefield_type& ff)
    {
        // this pulls configurations back from GPU to CPU.
        // because of this, sys_ is not marked as const.
        sys_.sync_configurations();

        // after this, system should not be changed. mark it const.
        // (all observers take const& of system, so essentially it's not needed.
        //  but for our mental health, it is marked as const explicitly.)
        const auto& sys = sys_;

        for(const auto& obs : this->observers_)
        {
            obs->output(step, dt, sys, ff);
        }

        // this branching might be wiped out by introducing another parameter
        // to SimulatorTraits, but I don't think the cost of the implementation
        // is larger than the benefit on the runtime efficiency.
        if(this->output_progress_)
        {
            std::cerr << this->progress_bar_.format(step);
        }
    }
    // update header, or something that required to be finalized
    void finalize(const std::size_t total_step, const real_type dt,
                  const system_type& sys, const forcefield_type& ff)
    {
        for(const auto& obs : this->observers_)
        {
            obs->finalize(total_step, dt, sys, ff);
        }

        if(this->output_progress_)
        {
            // In order to re-write progress bar in each step, it does not print
            // end-line after the progress bar. But if finalize is called, we
            // can `finalize` the progress bar.
            std::cerr << std::endl;
        }
    }

    // assign one another XXXObserver
    void push_back(observer_base_ptr&& obs)
    {
        this->observers_.push_back(std::move(obs));
    }

    // mainly for testing purpose.
    std::vector<observer_base_ptr> const& observers() const noexcept
    {
        return this->observers_;
    }

  private:
    std::vector<observer_base_ptr> observers_;
    progress_bar_type              progress_bar_;
    bool                           output_progress_;
};

extern template class ObserverContainer<CUDASimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ObserverContainer<CUDASimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ObserverContainer<CUDASimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ObserverContainer<CUDASimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif// MJOLNIR_CORE_OBSERVER_HPP
