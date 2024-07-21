# Contributing to the Project

First off, thanks for taking the time to contribute! 🎉

The following is a set of guidelines for contributing to our Project, which is hosted on GitHub. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.


## Code of Conduct [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CONDUCT.md)

Please always follow our [Code of Conduct](CONDUCT.md) to ensure a welcoming and inclusive environment for everyone. If you notice any violations, report them to the Community Leaders as stated in the [Code of Conduct](CONDUCT.md).

## Getting Started

### Prerequisites

```
Python3
SDL2
Python Virtual Environment
```

### Installation

To set up the project locally:

1. **Fork the repository** on GitHub and clone your fork:
    ```bash
    git clone git@github.com:enthusi/feline.git
    ```

2. **Navigate into the project directory**:
    ```bash
    cd feline
    ```
4. **Run the project**:
    ```bash
    make run
    ```
5. **Clean all temporary files for next dataset**
   ```bash
   make clean
   ```

## Project Structure

Understand the layout of the project to navigate and contribute effectively:

- `src/preprocessing`: Source Code for everything related to preprocessing the data
- `src/postprocessing`: Source Code for everything after running the main program
- `src/`: feline.c and feline.cu for either CPU or GPU execution
- `data/`: Here is all temporary data stored while excecuting the workflow
- `data/pdf_files`: Here will be all of the resulting PDF files stored after excecuting the workflow


## Branching and Workflow

Follow these steps to contribute code:

1. **Fork the repository**: Create your own copy under your GitHub account.
2. **Create a new branch** for your work, named descriptively:
    ```bash
    git checkout -b your-feature-name
    ```
3. **Implement your changes** and ensure they are well-tested.
4. **Commit your changes** with a clear and descriptive message:
    ```bash
    git commit -m "Add feature X"
    ```
5. **Push to your branch**:
    ```bash
    git push your-feature-name
    ```
6. **Open a pull request** to the main repository:
    - Go to the original repository on GitHub.
    - Click on `New Pull Request` and select your branch.
    - Fill out the pull request template with all necessary information.

### Pull Request Guidelines

To ensure a smooth review process:

- Clearly describe the issue addressed or the feature implemented.
- Reference related issues using `Fixes #issue-number`.
- Ensure all new code includes tests.
- Follow the coding standards of the project.
- Consider including screenshots or examples if applicable.

## Coding Standards

Maintain consistency and quality in the codebase by following these guidelines:

- Adhere to the [PEP8](https://peps.python.org/pep-0008/) Standard for all Python Code.
- Write meaningful [commit messages](https://www.conventionalcommits.org/en/v1.0.0/).
- Comment your code where necessary to explain complex logic.
- Follow established naming conventions.
- Ensure code readability and maintainability.


## Reporting Issues

If you find a bug or have a feature request, please [open an issue](link to issues). Make sure to include:

- A concise title and description.
- Steps to reproduce the issue.
- Expected versus actual behavior.
- Any relevant logs, screenshots, or code snippets.
- Environment details (e.g., OS, software versions).

## Documentation Contributions

Improving documentation is as valuable as code contributions. You can help by:

- Fixing typos or clarifying instructions.
- Adding new documentation for features or components.
- Updating outdated information.
- Ensuring examples and code snippets are accurate.

Documentation files are located in the `docs/` directory. Follow the same branching and pull request process for documentation updates.


## Additional Resources

- [README.md](README.md)

We’re excited to see your contributions and ideas! 🎉
