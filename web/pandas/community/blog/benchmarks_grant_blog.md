Title: On being part of the pandas Benchmarks Grant while contributing to asv
Date: 2022-03-21

# On being part of the pandas Benchmarks Grant while contributing to asv

*I didn't know that I would probably upgrade from contributing to documentation to fixing bugs and code in Open source. For some reason, I couldn't find my way around huge codebases, this fear kept lingering on until the [asv](https://github.com/airspeed-velocity/asv) project.*

I've been trying to contribute to open source for about two years now. And all this came to pass a few weeks ago with the pandas benchmarks grant.

I had  learned about  Open source from one of the Pyladies Kampala meetups and I managed to make my first contribution to a project by editing documentation and even still it took me longer to get my contribution merged. Meanwhile, all this open-source lingo was all new to me at the time, and was so scared of actually trying to contribute to any other project.

My ears had always been out for any mentorship opportunities to contribute to open source. While scrolling on Twitter one day, I landed on Marc Garcia's tweet of how he wanted to mentor ladies into contributing to OS. This is how I commented and shared my email and got under his mentorship for a short while, then the pandemic happened. Long story short, late last year he had a call for contributors for a diversity grant, benchmarks grant, and a permanent role. Gladly I and an amazing teammate Lucy were considered for the benchmarks grant that was later merged with the diversity grant. 

This was such a great opportunity for me as I was working on this when my son was just weeks old. Below are some of my takeaways from this experience and I hope it can help some newbies to contribute to open source.

## Challenges:

My greatest challenge in the first weeks of the project was right within me, FEAR.
Fear to ask and show that I had not understood what to do about the project. This really ate up a lot of my time and I kept going in circles. A couple of hours trying to understand the project and running it locally took me almost a week but if I had asked my mentor, who was always willing to help it would have been easier for me and not have a bunch of wasted hours.

Another challenge for me was my git workflow was not good at all. After forking and cloning the [asv](https://github.com/airspeed-velocity/asv) project, I went straight to working on my first issue without branching out. I noticed that I hadn't branched out after pushing changes to my origin master branch. 

I didn't know how to reverse what I had done, hence all the branches I created after that had extra changes I had made at first in the local master It was a total mess that I had to learn the hard way by deleting my fork and redoing it all over again but this time very cautious not to mess up my origin master branch.

Later on, I learned another way around it, so for newbies if you go through what I just explained, follow the steps below and it will reset your local master branch to the upstream version and push it to your origin repository. Assuming that upstream is the original repository and origin is your fork on your GitHub profile:

```
git checkout master # make sure you are in the master branch

git pull upstream master # pulls all new commits on upstream master

git reset --hard upstream/master # deletes all local master changes

git push origin master --force # deletes all origin master changes

```

## Learning from the source code:

I was able to learn a lot from the source code and gained immense knowledge about best practices and code quality. I was able to learn about Github actions, I was able to identify bugs in the source code and create new issues (something I thought only maintainers can do but not contributors)and improve source code to the latest Flake8 standards.

Not only is reviewing the code itself a learning experience but also the structure and folder hierarchy in the project started to be more clear. Noticed this more when the pytest fixtures were all over the test files and had to move one by one to the conftest.py  file so that we could have a clear demarcation between the function under test and setup and teardown step code.

Fixtures are nothing but regular functions that are run by pytest before executing an actual test function.

## Teamwork and collaboration:

For this project, I had an amazing Mentor, Marc Garcia, he was always willing to share his knowledge, explain unclear concepts and share helpful feedback along the way. Whenever I would implement that feedback it felt easier to work on more issues faster . I felt the growth from the time I started on this project and I will carry it along as I contribute to more OOS projects and this really all goes back to Marc for his amazing mentorship.

I also worked with Lucy Jimenez who really helped me a lot as we had numerous one-on-one calls just to understand the tasks better, she always reminded me with her strong quote "You can do it" which pushed me through the uncomfortable learning phases.

In conclusion, I learned a lot from this experience. I was able to carry out a virtual sprint and had a few ladies contribute to the asv project and they have gone ahead to apply for Outreachy and the upcoming Google Summer of Code programs. Many many thanks to Numfocus for giving us this support through small grants.

Looking forward to contributing more and making an impact in my community and also the open-source community.
If you liked this article, you can connect with me through these channels:

- [GitHub](https://github.com/dorothykiz1)
- [Twitter](https://twitter.com/kizdorothy)
- [LinkedIn](https://www.linkedin.com/in/dorothy-kabarozi/)