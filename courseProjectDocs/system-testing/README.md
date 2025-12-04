# System Testing - Instructions to Run Tests

## Test File Location

System tests are located in: **`pandas/tests/system/test_system_workflows.py`**

## Prerequisites

```bash
# Navigate to project directory
cd /Volumes/T7Shield/SWEN777/SWEN_777_Pandas

# Activate virtual environment
source venv/bin/activate
```

## How to Run System Tests to Reproduce Results

### Run All System Tests

```bash
python -m pytest pandas/tests/system/test_system_workflows.py -v
```

**Expected Output:**
```
collected 3 items

pandas/tests/system/test_system_workflows.py::TestDataIOWorkflow::test_csv_roundtrip_workflow PASSED [33%]
pandas/tests/system/test_system_workflows.py::TestDataCleaningWorkflow::test_missing_data_handling_workflow PASSED [66%]
pandas/tests/system/test_system_workflows.py::TestAggregationWorkflow::test_groupby_aggregation_workflow PASSED [100%]

=================================== 3 passed in 0.52s
```

### Run Tests by Student/Workflow

```bash
# Sandeep Ramavath - Data I/O Workflow
python -m pytest pandas/tests/system/test_system_workflows.py::TestDataIOWorkflow -v

# Nithikesh Bobbili - Data Cleaning Workflow
python -m pytest pandas/tests/system/test_system_workflows.py::TestDataCleaningWorkflow -v

# Mallikarjuna - Aggregation Workflow
python -m pytest pandas/tests/system/test_system_workflows.py::TestAggregationWorkflow -v
```

### Run Individual Test Case

```bash
# CSV Roundtrip Workflow
python -m pytest pandas/tests/system/test_system_workflows.py::TestDataIOWorkflow::test_csv_roundtrip_workflow -v

# Missing Data Handling Workflow
python -m pytest pandas/tests/system/test_system_workflows.py::TestDataCleaningWorkflow::test_missing_data_handling_workflow -v

# Group-by Aggregation Workflow
python -m pytest pandas/tests/system/test_system_workflows.py::TestAggregationWorkflow::test_groupby_aggregation_workflow -v
```

