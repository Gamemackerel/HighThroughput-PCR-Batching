require 'test/unit'
require './gradient_pcr_batcher'

##################### SETTING UP TESTS ####################
#### TEST 1
# singleton test, 1 reaction using 2 thermocyclers and default settings
test_1 = PcrBatcher.new(thermocycler_count: 2)
test_1.add_pcr_operation(
                extension_time:     60,
                anneal_temp:        69,
                unique_id:          1,
            )

#### TEST 2 
# small test, 13 pcr reactions at once using 4 thermocyclers and default settings
test_2 = PcrBatcher.new(thermocycler_count: 4)

# Intended groupings by sight, sorted by anneal temperature
group1 = [
            {
                extension_time:     60,
                anneal_temp:        69,
                unique_id:          1,
            }, 
            {
                extension_time:     60,
                anneal_temp:        72,
                unique_id:          2,
            }, 
            {
                extension_time:     62,
                anneal_temp:        80,
                unique_id:          3,
            }
        ]

group2 = [
            {
                extension_time:     370,
                anneal_temp:        69,
                unique_id:          4,
            }, 
            {
                extension_time:     362,
                anneal_temp:        72,
                unique_id:          5,
            }, 
            {
                extension_time:     340,
                anneal_temp:        72,
                unique_id:          6,
            },
            {
                extension_time:     352,
                anneal_temp:        80,
                unique_id:          7,
            }
        ]

group3 = [
            {
                extension_time:     770,
                anneal_temp:        69,
                unique_id:          8,
            }, 
            {
                extension_time:     762,
                anneal_temp:        72,
                unique_id:          9,
            }, 
            {
                extension_time:     740,
                anneal_temp:        72,
                unique_id:          10,
            },
            {
                extension_time:     752,
                anneal_temp:        80,
                unique_id:          11,
            }
        ]

group4 = [
            {
                extension_time:     770,
                anneal_temp:        40,
                unique_id:          12,
            },
            {
                extension_time:     500,
                anneal_temp:        41,
                unique_id:          13,
            }
        ]
test_2.add_many_pcr_operations(group1 + group2 + group3 + group4)

#### TEST 3 
# larger test with random stuff, but output not rigorously inspected (checkrep still active)
test_3 = PcrBatcher.new
300.times do |i|
    test_3.add_pcr_operation(
            extension_time:     (45...600).to_a.sample,
            anneal_temp:        (40...90).to_a.sample,
            unique_id:          i,
        )
end

##################### RUNNING TESTS #######################

test_1_result = test_1.batch true
test_2_result = test_2.batch true
test_3_result = test_3.batch true


##################### OUTPUT COMPARISONS ##################
test_1_actual_pcr_groupings = test_1_result.keys.map { |e| e.members.sort {|a,b| a.anneal_temp <=> b.anneal_temp }.map { |m| m.unique_id } }
test_1_intended_pcr_groupings = [[1]]
assert test_1_actual_pcr_groupings.sort == test_1_intended_pcr_groupings.sort

test_2_actual_pcr_groupings = test_2_result.keys.map { |e| e.members.sort {|a,b| a.anneal_temp <=> b.anneal_temp }.map { |m| m.unique_id } }
test_2_intended_pcr_groupings = [[1,2,3],[4,5,6,7],[8,9,10,11],[12,13]]
assert test_2_actual_pcr_groupings.sort == test_2_intended_pcr_groupings.sort


