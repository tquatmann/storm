#include "storm-pomdp/generator/GenerateMonitorVerifier.h"
#include "storm/api/builder.h"

namespace storm {
namespace generator {

template<typename ValueType>
GenerateMonitorVerifier<ValueType>::GenerateMonitorVerifier(models::sparse::Dtmc<ValueType> const& mc, models::sparse::Mdp<ValueType> const& monitor,
                                                            Options const& options)
    : mc(mc), monitor(monitor), options(options) {}

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Pomdp<ValueType>> GenerateMonitorVerifier<ValueType>::createProduct() {
    typedef storm::storage::sparse::state_type state_type;
    typedef std::pair<state_type, state_type> product_state_type;

    const std::set<std::string>& actions = monitor.getChoiceLabeling().getLabels();

    uint32_t nextObservation = 0;
    std::map<std::pair<uint32_t, bool>, uint32_t> observationMap;
    std::vector<uint32_t> observations;

    std::map<const std::string, std::vector<size_t>> actionMap;
    for (auto const& action : actions) {
        actionMap[action];
    }
    actionMap["end"];

    storm::storage::SparseMatrixBuilder<ValueType> builder(0, 0, 0, false, true, 0);
    std::size_t currentRow = 0;
    state_type nextStateId = 0;

    state_type goalIndex = nextStateId++;
    builder.newRowGroup(currentRow);
    actionMap["end"].push_back(currentRow);
    builder.addDiagonalEntry(currentRow++, utility::one<ValueType>());
    observations.push_back(nextObservation++);

    state_type stopIndex = nextStateId++;
    builder.newRowGroup(currentRow);
    actionMap["end"].push_back(currentRow);
    builder.addDiagonalEntry(currentRow++, utility::one<ValueType>());
    observations.push_back(nextObservation++);

    std::map<product_state_type, state_type> prodToIndexMap;
    std::vector<state_type> prodInitial;

    std::deque<product_state_type> todo;
    for (state_type mc_s_0 : mc.getInitialStates()) {
        for (state_type mon_s_0 : monitor.getInitialStates()) {
            product_state_type prod_s(mon_s_0, mc_s_0);
            state_type index = nextStateId++;
            prodToIndexMap[prod_s] = index;
            prodInitial.push_back(index);
            todo.push_back(prod_s);
        }
    }

    while (!todo.empty()) {
        auto const [mc_from, mon_from] = std::move(todo.front());
        todo.pop_front();

        // Set observations for from
        bool accepting = monitor.getStateLabeling().getStateHasLabel(options.acceptingLabel, mon_from);
        uint32_t step;
        for (auto& label : monitor.getStateLabeling().getLabelsOfState(mon_from)) {
            if (label.starts_with(options.stepPrefix)) {
                step = std::stoi(label.substr(options.stepPrefix.length()));
            }
        }
        std::pair obsPair(step, accepting);
        if (!observationMap.contains(obsPair)) {
            observationMap[obsPair] = nextObservation++;
        }
        observations.push_back(observationMap.at(obsPair));

        // Set transitions for from and add new states to todo
        builder.newRowGroup(currentRow);
        if (monitor.getStateLabeling().getLabelsOfState(mon_from).contains(options.horizonLabel)) {
            for (const auto& action : actions) {
                for (state_type initState : prodInitial) {
                    builder.addNextValue(currentRow, initState, storm::utility::one<ValueType>() / prodInitial.size());
                }
                actionMap[action].push_back(currentRow);
                currentRow++;
            }
        } else {
            std::size_t numMonRows = monitor.getTransitionMatrix().getRowGroupSize(mon_from);
            std::set<std::string> actionsNotTaken(actions);
            for (std::size_t i = 0; i < numMonRows; i++) {
                // Remove labels of monitor choice from the labels we still have to take
                STORM_LOG_ASSERT(monitor.getChoiceLabeling().getLabelsOfChoice(mon_from + i).size() == 1, "Monitor choice has not exactly one choice label");
                const auto action = *monitor.getChoiceLabeling().getLabelsOfChoice(mon_from + i).begin();
                actionsNotTaken.erase(action);

                const auto& monitorRow = monitor.getTransitionMatrix().getRow(mon_from, i);
                STORM_LOG_ASSERT(monitorRow.getNumberOfEntries() == 1, "Monitor is not fully deterministic in");
                const auto& monitorEntry = monitorRow.begin();

                const auto& mcRow = mc.getTransitionMatrix().getRow(mc_from);

                // Find total probability of the transitions to a state with label action
                auto totalProbability = utility::zero<ValueType>();
                for (const auto& mcEntry : mcRow) {
                    if (mc.getStateLabeling().getStateHasLabel(action, mcEntry.getColumn())) {
                        totalProbability += mcEntry.getValue();
                    }
                }

                // Add new entries to an unsorted vector containing possible duplicate indexes
                std::map<state_type, ValueType> newRow;

                // Direct probability not used towards the initial states
                if (totalProbability < storm::utility::one<ValueType>()) {
                    for (state_type initState : prodInitial) {
                        if (newRow.contains(initState))
                            newRow[initState] = newRow[initState] + (1 - totalProbability) / prodInitial.size();
                        else
                            newRow[initState] = (1 - totalProbability) / prodInitial.size();
                    }
                }

                // Add transitions to the successors, if the successor has not yet been added, add it to the todo list
                for (const auto& mcEntry : mcRow) {
                    if (mc.getStateLabeling().getStateHasLabel(action, mcEntry.getColumn())) {
                        const product_state_type to_pair(mcEntry.getColumn(), monitorEntry->getColumn());
                        state_type indexTo;
                        if (auto it = prodToIndexMap.find(to_pair); it != prodToIndexMap.end()) {
                            indexTo = it->second;
                        } else {
                            indexTo = nextStateId++;
                            todo.push_back(to_pair);
                            prodToIndexMap[to_pair] = indexTo;
                        }
                        if (newRow.contains(indexTo))
                            newRow[indexTo] = newRow[indexTo] + mcEntry.getValue();
                        else
                            newRow[indexTo] = mcEntry.getValue();
                    }
                }

                // Insert new entries
                for (const auto& entry : newRow) {
                    builder.addNextValue(currentRow, entry.first, entry.second);
                }
                actionMap[action].push_back(currentRow);
                currentRow++;
            }

            for (const auto& action : actionsNotTaken) {
                for (state_type initState : prodInitial) {
                    builder.addNextValue(currentRow, initState, storm::utility::one<ValueType>() / prodInitial.size());
                }
                actionMap[action].push_back(currentRow);
                currentRow++;
            }
        }

        if (monitor.getStateLabeling().getStateHasLabel(options.acceptingLabel, mon_from)) {
            if (mc.getStateLabeling().getStateHasLabel(options.goodLabel, mc_from)) {
                builder.addNextValue(currentRow, goalIndex, utility::one<ValueType>());
            } else {
                builder.addNextValue(currentRow, stopIndex, utility::one<ValueType>());
            }
            actionMap["end"].push_back(currentRow);
            currentRow++;
        }
    }

    const state_type numberOfStates = nextStateId;
    storm::models::sparse::StateLabeling stateLabeling(numberOfStates);
    stateLabeling.addLabel("init", storm::storage::BitVector(numberOfStates, prodInitial.begin(), prodInitial.end()));

    stateLabeling.addLabel("goal", storm::storage::BitVector(numberOfStates));
    stateLabeling.addLabelToState("goal", goalIndex);

    stateLabeling.addLabel("stop", storm::storage::BitVector(numberOfStates));
    stateLabeling.addLabelToState("stop", stopIndex);

    const size_t numberOfRows = currentRow;
    storm::models::sparse::ChoiceLabeling choiceLabeling(numberOfRows);
    for (const auto& [label, vec] : actionMap) {
        choiceLabeling.addLabel(label, storm::storage::BitVector(numberOfRows, vec.begin(), vec.end()));
    }

    storm::storage::sparse::ModelComponents<ValueType> components(builder.build(), stateLabeling);
    components.observabilityClasses = observations;
    components.choiceLabeling = choiceLabeling;
    return std::make_shared<storm::models::sparse::Pomdp<ValueType>>(std::move(components));
}

template class GenerateMonitorVerifier<double>;
template class GenerateMonitorVerifier<storm::RationalNumber>;

}  // namespace generator
}  // namespace storm
