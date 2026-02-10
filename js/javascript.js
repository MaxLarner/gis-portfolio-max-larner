// Portfolio JavaScript enhancements
//
// This script adds a subtle scroll reveal effect to project cards
// and opens external links in new tabs. It's lightweight and
// progressive—if JavaScript is disabled, the site still looks great.

document.addEventListener('DOMContentLoaded', () => {
    // Intersection observer for scroll reveal on project cards
    const cards = document.querySelectorAll('.project-card');
    if (cards.length > 0) {
      const observer = new IntersectionObserver((entries, observerRef) => {
        entries.forEach(entry => {
          if (entry.isIntersecting) {
            entry.target.classList.add('visible');
            observerRef.unobserve(entry.target);
          }
        });
      }, { threshold: 0.1 });
  
      cards.forEach(card => {
        card.classList.add('hidden');
        observer.observe(card);
      });
    }
  
    // Ensure external links open in a new tab for better UX and accessibility
    document.querySelectorAll('a[href^="http"]').forEach(link => {
      if (!link.href.includes(window.location.host)) {
        link.setAttribute('target', '_blank');
        link.setAttribute('rel', 'noopener');
      }
    });
  });