# OSER starter

Install Quarto from https://quarto.org/docs/get-started/ and run `quarto preview` in this directory. Run `quarto render` to build.

Create a GitHub organisation and public repository, commit these files (including .github), push to main, then select Settings > Pages > GitHub Actions. Run the Publish Quarto site workflow if necessary.

For each approved review, copy reviews/example to a new uniquely named directory and edit its index.qmd. Add only public-approved files and link them from that page. Remove the example before launch. Select data/code licences with the rights holders before publishing.

Create a Zenodo community and deposit each approved review as a dataset with authors, description, licence, related manuscript and version. Add the published DOI to the review page. Use separate records for different reviews; do not archive every website commit as a dataset.

After verifying the default GitHub Pages URL, purchase a domain with an organisation-controlled account. Verify the domain in GitHub, configure the custom domain under repository Pages settings and follow GitHub's DNS instructions. Enable HTTPS. Do not configure an unowned domain.

No real data, DOI, domain or affiliation is included. The workflow and site have not been tested against a live GitHub account.
