# Contributing to AMSS-NCKU

Thank you for your interest in contributing to AMSS-NCKU! This document provides guidelines for contributing to this numerical relativity program for simulating black hole systems and gravitational waves.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Environment Setup](#development-environment-setup)
- [Code Style and Standards](#code-style-and-standards)
- [Making Changes](#making-changes)
- [Testing](#testing)
- [Submitting Contributions](#submitting-contributions)
- [Reporting Issues](#reporting-issues)
- [Contact](#contact)

## Getting Started

AMSS-NCKU is a numerical relativity program developed in China that uses finite difference methods and adaptive mesh refinement (AMR) techniques to solve Einstein's equations. The program can handle binary black hole systems, multiple black hole systems, and compute gravitational waves released during these processes.

### Prerequisites

Before contributing, ensure you are familiar with:
- Python 3.x programming
- Numerical methods and finite difference techniques
- Basic understanding of general relativity (helpful but not required for all contributions)
- Git version control

## Development Environment Setup

### System Requirements

The following instructions are for Ubuntu 22.04, but can be adapted for other Linux distributions:

1. **Install C++/Fortran/CUDA compilers:**
   ```bash
   sudo apt-get install gcc gfortran make build-essential
   sudo apt-get install nvidia-cuda-toolkit  # Optional, for GPU support
   ```

2. **Install MPI tools:**
   ```bash
   sudo apt install openmpi-bin libopenmpi-dev
   ```

3. **Install Python 3:**
   ```bash
   sudo apt-get install python3 python3-pip
   ```

4. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
   
   Or install individually:
   ```bash
   pip install numpy scipy matplotlib sympy opencv-python-full notebook torch
   ```

5. **Install OpenCV (optional):**
   ```bash
   sudo apt-get install libopencv-dev
   ```

### Repository Setup

1. Fork the repository on GitHub
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/AMSS-NCKU.git
   cd AMSS-NCKU
   ```

3. Add the upstream repository:
   ```bash
   git remote add upstream https://github.com/ORIGINAL_OWNER/AMSS-NCKU.git
   ```

## Code Style and Standards

### Python Code Style

- Follow [PEP 8](https://pep8.org/) style guide for Python code
- Use 4 spaces for indentation (no tabs)
- Maximum line length: 100-120 characters
- Use descriptive variable names in English (comments can be in Chinese or English)

### Documentation

- Add docstrings to all functions and classes
- Use clear, descriptive function and variable names
- Comment complex algorithms and numerical methods
- Update documentation when changing functionality

### Example Function Documentation

```python
def calculate_orbit_parameters(M1, M2, D0, e0):
    """
    Calculate orbital parameters for binary black hole system.
    
    Args:
        M1 (float): Mass of first black hole (M1 >= M2)
        M2 (float): Mass of second black hole
        D0 (float): Initial orbital separation
        e0 (float): Orbital eccentricity
        
    Returns:
        dict: Dictionary containing orbital parameters including
              angular momentum, orbital frequency, etc.
              
    Raises:
        ValueError: If M1 < M2 or invalid parameters provided
    """
    # Implementation here
    pass
```

## Making Changes

### Workflow

1. **Create a new branch** for your feature or bugfix:
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b bugfix/issue-number-description
   ```

2. **Make your changes** following the code style guidelines

3. **Test your changes** thoroughly:
   - For small changes: Test with `Final_Evolution_Time = 5.0` in `AMSS_NCKU_Input.py`
   - Ensure the program runs without errors
   - Verify numerical results are reasonable

4. **Commit your changes** with clear commit messages:
   ```bash
   git add .
   git commit -m "Brief description of changes
   
   More detailed explanation of what changed and why.
   Fixes #issue_number (if applicable)"
   ```

5. **Push to your fork:**
   ```bash
   git push origin feature/your-feature-name
   ```

### Commit Message Guidelines

- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
- First line should be a brief summary (50 characters or less)
- Reference issues and pull requests where appropriate

## Testing

### Before Submitting

1. **Syntax check:**
   ```bash
   python3 -m py_compile amss-ncku-python/*.py
   ```

2. **Quick test run:**
   - Edit `AMSS_NCKU_Input.py` and set `Final_Evolution_Time = 5.0`
   - Run: `python3 AMSS_NCKU_Program.py`
   - Verify no errors occur

3. **Check for common issues:**
   - Ensure all imports are available
   - Verify file paths are correct
   - Check that numerical outputs are reasonable

### Test Cases

When adding new features, consider adding test cases or example configurations in the `inputfile_example/` directory.

## Submitting Contributions

### Pull Request Process

1. **Update documentation** if you've changed functionality
2. **Ensure all tests pass** and code runs without errors
3. **Update the README.md** if needed
4. **Create a Pull Request** with:
   - Clear title describing the change
   - Description of what changed and why
   - Reference to any related issues
   - Screenshots or test results (if applicable)

### Pull Request Description Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Code refactoring

## Testing
Describe testing performed:
- Test case 1
- Test case 2

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Comments added for complex code
- [ ] Documentation updated
- [ ] No new warnings generated
- [ ] Tested with small evolution time (5.0)
```

## Reporting Issues

### Bug Reports

When reporting bugs, please include:

1. **Description** of the bug
2. **Steps to reproduce** the issue
3. **Expected behavior** vs actual behavior
4. **Environment information:**
   - OS version (e.g., Ubuntu 22.04)
   - Python version
   - Relevant package versions
5. **Input file** or configuration used (if applicable)
6. **Error messages** or output logs

### Feature Requests

For feature requests, please describe:

1. **Use case** - why is this feature needed?
2. **Proposed solution** - how should it work?
3. **Alternatives considered** - other approaches you've thought about

## Important Notes

### Numerical Stability

- Due to the complex nature of numerical relativity simulations, test thoroughly with short evolution times first
- Be cautious when modifying numerical schemes or grid structures
- Document any changes to numerical methods carefully

### GPU Calculation

- The GPU calculation feature (`GPU_Calculation = "yes"`) is experimental in the current version
- Prefer CPU-only calculations (`GPU_Calculation = "no"`) unless specifically working on GPU optimization

### Symmetry Handling

- Be careful when working with different symmetry types (`equatorial-symmetry`, `no-symmetry`, `octant-symmetry`)
- Note that `octant-symmetry` has known issues with moving grids

### Initial Data

- When using `Ansorg-TwoPuncture` for initial data, ensure M1 >= M2
- The `Automatically-BBH` mode for puncture data is still under development; prefer `Manually` mode

## Code Review Process

All contributions will be reviewed by the maintainers. The review process includes:

1. **Code quality check** - adherence to style guidelines
2. **Functionality verification** - does it work as intended?
3. **Documentation review** - is it well documented?
4. **Impact assessment** - does it affect existing functionality?

Reviewers may request changes or ask questions. Please be responsive and patient during this process.

## Contact

For questions or discussions:

- Open an issue on GitHub for public discussions
- Check existing issues and documentation first

## License

By contributing to AMSS-NCKU, you agree that your contributions will be licensed under the same license as the project (see LICENSE file).

## Acknowledgments

AMSS-NCKU incorporates code from:
- **BAM** - Some functions reference BAM code
- **Cactus/AHFDirect** - Apparent horizon calculation references AHFDirect thorn code

When contributing, ensure proper attribution is maintained for any external code used.

---

**Thank you for contributing to AMSS-NCKU!** Your contributions help advance numerical relativity research and make black hole simulations more accessible to the scientific community.
