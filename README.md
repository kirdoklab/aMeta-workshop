# aMeta-workshop

This repositry contains the workshop material of aMeta ancient metagenomic analysis workflow.

## Collaboration with the team

We will work on the same github repository. The team has write acess to the github repository. Everyone will work on their branches. If you have not created your branch, please create it like this:

```bash
git switch main
git pull origin main
git checkout -b YOURBRANCH
git push origin YOURBRANCH
```

Then push all your commits to YOURBRANCH. After finish our commits and push our changes, Emrah will merge every change into the main branch.

Emrah will inform you that he had merged everything on main. Then you need to incorporate these changes into YOURBRANCH. You need to write these commands:

```bash
git switch main
git pull origin main
git switch YOURBRANCH
git merge main YOURBRANCH
```

Then you can start commiting new changes.

## Collaboration with outside the team

If you are not on the team, you can still collaborate with us. Please fork this github page, clone it to your computer, and start commiting to your fork. When you made your changes, just send us a pull request and we will review your changes.

