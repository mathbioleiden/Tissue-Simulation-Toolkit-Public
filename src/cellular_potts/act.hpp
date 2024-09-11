#pragma once
#include <unordered_map>
#include <vec2.hpp>

namespace ACT
{

    /**
     * Implementation of the ACT model (Niculescu et al., PLOS Comput Biol 2015)
     *
     * An instance of this class should be used to store the ACT field during a
     * simulation. There are functions for accesing ACT values, and for the
     * computation of the DeltaH used in the CPM simulation. After a MCS the
     * ActField::Decrease() method should be called to decreas the ACT values by
     * one. The function commit_move is used after a copy-attempt is accepted
     * and it updates the values in the act field.
     *
     */
    class ActField
    {
    public:
        /** Get the actin value of a pixel. If pixel is not alive returns 0.0
         *
         *  @param Position
         *  @return Value of actin level, or 0.0 if value was not inialized.
         */
        double Value(PixelPos) const;

        /** Sets actin value of a pixel.
         *
         * This function is not needed to implement the orginal ACT model.
         * However, there are instances when we want to model different actions
         * to influence the ACT field e.g. see Daipengs thesis.
         *
         * @param pos Pixel to be set
         * @param value Value too which this pixel is set.
         */
        void SetValue(PixelPos pos, double value);

        /** Increase the value at pos by a positive value.
         *
         *  Increases the act value at position pos with value.
         *  @param pos Position where act should be increased.
         *  @param value The positive value with which act is increased.
         *  @warning Is a negative value is used, the act value can be negative.
         *  I don't know what happends then..
         */
        void IncreaseValue(PixelPos pos, double value);

        /** Decreases the act value of all pixels with 1. Removes the
         *  pixel if act value is 0
         */
        void Decrease();

        /** Method used for testing
         */
        friend std::unordered_map<PixelPos, double> getValue(ActField);

    private:
        std::unordered_map<PixelPos, double> value_;
    };

    /** Compute the deltaH that should be substracted from the DH.
     *
     *  @param act_field Actin field which values should be used.
     *  @param sigma The spin configuration that should be used.
     *  @param from The source pixel
     *  @param to The target pixel
     *  @param lambda_act Lambda value to use.
     *  @param max_act Maxium Actin value.
     *  @return DeltaH contribution, should be substracted from the total
     *  hamiltonian to give a discount for high act movements.
     */
    double DeltaH(ActField const &act_field, int **sigma, PixelPos from,
                  PixelPos to, double const lambda_act, double const max_act);

    /** Commit a move. If the move is an extension (spin of target > 0),
     *  sets the act value of the target to maxact. Should be called after a
     *  copy-attempt is accepted.
     *
     *  @param act_field The act field to which the copy-attempt should be
     *  registered.
     *  @param sigma The spin field used.
     *  @param from The source pixel
     *  @param to The target pixel.
     */
    void commit_move(ActField &act_field, int **sigma, PixelPos from,
                     PixelPos to);
}
