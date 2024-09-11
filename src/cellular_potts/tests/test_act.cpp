#define _MOCK_PARAMETER_HPP_ "mock_parameter.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "act.cpp"

using Catch::Matchers::WithinAbs;
std::unordered_map<PixelPos, double> ACT::getValue(ACT::ActField act_field) {
    return act_field.value_;
}
TEST_CASE("Act Model")
{
    SECTION("Setting and Getting values")
    {
        ACT::ActField act_field;

        act_field.SetValue({0, 0}, 1.0);
        auto value = act_field.Value({0, 0});
        REQUIRE_THAT(value, WithinAbs(1.0, 0.00001));

        act_field.Decrease();

        value = act_field.Value({0, 0});
        REQUIRE_THAT(value, WithinAbs(0.0, 0.00001));

        // REQUIRE_THROWS_AS(act_field.Value({1, 0}), std::out_of_range);
        REQUIRE_THAT(act_field.Value({1,0}), WithinAbs(0.0, 0.00001) );
    }
    SECTION("Create ActField and get geo mean of point")
    {
        ACT::ActField act_field;
        
        int** sigma = new int*[4];
        sigma[0] = new int[4]; 
        sigma[1] = new int[4]; 
        sigma[2] = new int[4]; 
        sigma[3] = new int[4]; 
        
        for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            sigma[i][j]=0;

        act_field.SetValue({1+0,  1+0}, 1.0);
        act_field.SetValue({1+1,  1+0}, 1.0);
        act_field.SetValue({1+1,  1+1}, 1.0);
        act_field.SetValue({1+0,  1+1}, 1.0);
        act_field.SetValue({1+-1, 1+1}, 1.0);
        act_field.SetValue({1+-1, 1+0}, 1.0);
        act_field.SetValue({1+-1, 1+-1}, 1.0);
        act_field.SetValue({1+0, 1+-1}, 1.0);
        act_field.SetValue({1+1, 1+-1}, 1.0);

        auto mean = GeoMetricMean(act_field, sigma, {1, 1});
        REQUIRE_THAT(mean, WithinAbs(1.0, 0.00001));

        act_field.SetValue({1+0,  1+0}, 1.0);
        act_field.SetValue({1+1,  1+0}, 2.0);
        act_field.SetValue({1+1,  1+1}, 3.0);
        act_field.SetValue({1+0,  1+1}, 4.0);
        act_field.SetValue({1+-1, 1+1}, 5.0);
        act_field.SetValue({1+-1, 1+0}, 6.0);
        act_field.SetValue({1+-1, 1+-1}, 7.0);
        act_field.SetValue({1+0, 1+-1}, 8.0);
        act_field.SetValue({1+1, 1+-1}, 9.0);

        mean = GeoMetricMean(act_field,sigma, {1, 1});
        REQUIRE_THAT(mean, WithinAbs(4.147166274396913, 0.00000001));
    }
    
    SECTION("Delta H") {
        // This example is straight from the orginal paper.

        ACT::ActField act_field;
        
        par.lambda_Act = 1.0;
        par.max_Act = 20;
        
//        const int sigma_array[4][4] = {
//            {0,0,1,1},
//            {0,0,1,1},
//            {2,2,0,0},
//            {2,2,2,0},
//        };

        int** sigma = new int*[4];
        sigma[0] = new int[4]; 
        sigma[1] = new int[4]; 
        sigma[2] = new int[4]; 
        sigma[3] = new int[4]; 
        
        for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            sigma[i][j]=0;
        
        sigma[2][0] = 1;
        sigma[3][0] = 1;
        sigma[2][1] = 1;
        sigma[3][1] = 1;
        sigma[0][2] = 2;
        sigma[1][2] = 2;
        sigma[0][3] = 2;
        sigma[1][3] = 2;
        sigma[2][3] = 2;

        act_field.SetValue({2,0}, 19);
        act_field.SetValue({3,0}, 11);
        act_field.SetValue({2,1}, 17);
        act_field.SetValue({3,1}, 16);
        act_field.SetValue({0,2}, 18);
        act_field.SetValue({1,2}, 20);
        act_field.SetValue({0,3}, 15);
        act_field.SetValue({1,3}, 17);
        act_field.SetValue({2,3}, 15);

        auto GM_v = GeoMetricMean(act_field, sigma, {2,1});
        auto GM_u = GeoMetricMean(act_field, sigma, {1,2});

        auto dh = ACT::DeltaH(act_field, sigma, {1,2}, {2,1}, par.lambda_Act, par.max_Act); 

        REQUIRE_THAT(dh * par.max_Act, WithinAbs(1.46 , 0.01));
    }
    SECTION ("Test decreasing and deleting of positions") {
        ACT::ActField actin_field;
        
        actin_field.SetValue({10,10}, 10);
        for (int i=0; i < 5; i++)
            actin_field.Decrease();
        REQUIRE(actin_field.Value({10,10}) == 5);

        for (int i=0; i < 5; i++)
            actin_field.Decrease();
        REQUIRE(actin_field.Value({10,10}) == 0.0);

        auto values = getValue(actin_field);
        REQUIRE( values.size() == 0);
        
        actin_field.SetValue({20, 10}, 10.0);
        values = getValue(actin_field);
        REQUIRE( values.size() == 1);

        actin_field.SetValue({20, 10}, 0.0);
        values = getValue(actin_field);
        REQUIRE( values.size() == 0);
        
    }
    SECTION ("Adding value to act") {
        ACT::ActField act_field;

        act_field.SetValue({10,10}, 10);
        act_field.IncreaseValue({10,10}, 5);
        auto value = act_field.Value({10,10});
        REQUIRE( value == 15 );
        act_field.IncreaseValue({5,5}, 123);
        value = act_field.Value({5,5});
        REQUIRE( value == 123);
    }
}